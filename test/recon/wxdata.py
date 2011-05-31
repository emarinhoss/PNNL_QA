r"""Provides an interface to read and WarpX hyperbolic solver data. The
numpy and PyTables package is required. The code works with 1D/2D/3D
simulations.

To read data from a run with input file 'axishock-wv.inp' use::

Reading a dump file
-------------------

To read data dump 1 with run prefix 'axishock-wv', for example,
use::

 >>> td = wxdata.WxData('axishock-wv', 5)

At this point no updatable data is actually read into memory but
meta-data about the simulation and the dump can be accessed. An
additional parameter can be passed which specifies the ``flavor`` of
the arrays. If specified this should be one of ``numpy`` (default) or
``numeric``.  The time at which the data was written can be found by::

 >>> td.time
 0.35

The top level combosolver can be determined::

 >>> td.comboSolver()
 'comboSolver'

Determining variables and getting data from arrays
--------------------------------------------------

To determine the variables defined in the simulation use::

 >>> td.variables('comboSolver')
 ['q', 'qnew']

 We can now fetch data from a variable by::

  >>> q = td.read('qnew')

The variable ``q`` is of type ``WxArray``. At this stage no data is
actually read from the HDF5 file. One can determine some properties of
the array::
 
 >>> q.numComponents
 5L

gives the number of components per cell in the array,::

 >>> q.shape
 (600L, 400L)

gives the shape of the array. The grid on which the array lives can
be accessed::

 >>> grid = q.grid

Its properties can be determined::

 >>> grid.lowerBounds
 [0.0, 0.0]

 >>> grid.upperBounds
 [1.5, 1.0]

 >>> grid.numPhysCells
 [600, 400]

 >>> grid.dx
 [0.0025, 0.0025]

To actually read the data into python::

 >>> q1 = q.load()

which will load the full array into memory. One also load only a slice
of the data::

 >>> q2 = q[:,400,0]

which will read all values along the X-axis in the mid-plain. Note
that while slicing only that part of the array is read into memory and
the whole data is *not* loaded into memory.

Reading  an array from DG simulation
------------------------------------

To read a data array from a DG simulation, specify the name of the
aray you want to read and the number of equations::

 >>> q = d.readDG('qnew', 5)

which will read in the DG data from the file and perform the
interpolation to the DG quadrature points. If the DG scheme has
spatial order 2 and was run on a grid of 50x50, then the shape of q is
100x100x5. The X coordinates of each quadrature point is in q.grid.X
and the Y coordinates of the quadrature points is q.grid.Y. These can
be used to plot the solution as needed.

"""

import warnings
warnings.simplefilter('ignore', UserWarning)

import tables
import numpy
import os
import math

class WxInterpGrid:
    def __init__(self, X, Y=None, Z=None):
        self.ndims = 1
        self.X = X
        self.Y = Y
        self.Z = Z
        if self.Y: ndims = 2
        if self.Z: ndims = 3

class WxRange:
    r"""
    Represents a range of indices.
    """

    def __init__(self, lst):
        self.ndims = len(lst)/2
        self.lower = self.ndims*[0]
        self.upper = self.ndims*[0]
        self.length = self.ndims*[0]        
        for i in range(self.ndims):
            self.lower[i] = lst[2*i]
            self.upper[i] = lst[2*i+1]
            self.length[i] = self.upper[i] - self.lower[i]

class WxColIndexer:
    r"""
    Provides a class to do column-major indexing for arbitrary
    dimensional arrays.
    """

    def __init__(self, r):
        self.ndims = r.ndims
        # allocate memory for coefficients
        self.ai = (self.ndims+1)*[0]

        # set a_1, ..., a_N
        self.ai[1] = 1
        for i in range(2,self.ndims+1):
            self.ai[i] = self.ai[i-1]*r.length[i-2]

        # set a_0
        sm = 0
        for i in range(1,self.ndims+1):
            sm = sm + self.ai[i]*r.lower[i-1]
        self.ai[0] = -sm

    def index1(self, k1):
        return self.ai[0] + k1

    def index2(self, k1, k2):
        return self.ai[0] + k1 + self.ai[2]*k2

    def index3(self, k1, k2, k3):
        return self.ai[0] + k1 + self.ai[2]*k2 + self.ai[3]*k3;

    def index4(self, k1, k2, k3, k4):
        return self.ai[0] + k1 + self.ai[2]*k2 + self.ai[3]*k3 + self.ai[4]*k4

class WxGrid:
    r"""WxGrid(lower, upper, cells) -> WxGrid

    Returns an object representing a rectangular grid in WarpX
    """
    
    def __init__(self, lowerBounds, upperBounds, numPhysCells):
        self.lowerBounds = list(lowerBounds)
        self.upperBounds = list(upperBounds)
        self.numPhysCells = list(numPhysCells)
        self.ndims = len(lowerBounds)
        # construct cell spacing
        self.dx = self.ndims*[0.0]
        for i in range(self.ndims):
            self.dx[i] = (upperBounds[i]-lowerBounds[i])/numPhysCells[i]

        # construct axes for plotting

        # X-axis
        self.X = numpy.linspace(
            lowerBounds[0], upperBounds[0]-self.dx[0], numPhysCells[0]) + 0.5*self.dx[0]

        if (self.ndims == 2):
            # Y-axis
            self.Y = numpy.linspace(
                lowerBounds[1], upperBounds[1]-self.dx[1], numPhysCells[1]) + 0.5*self.dx[1]

        if (self.ndims == 3):
            # Y-axis
            self.Y = numpy.linspace(
                lowerBounds[1], upperBounds[1]-self.dx[1], numPhysCells[1]) + 0.5*self.dx[1]
            # Z-axis
            self.Z = numpy.linspace(
                lowerBounds[2], upperBounds[2]-self.dx[2], numPhysCells[2]) + 0.5*self.dx[2]

    def __repr__(self):
        rs = ''
        rs = rs + ' lowerBounds = %s \n' % str(self.lowerBounds)
        rs = rs + ' upperBounds = %s \n' % str(self.upperBounds)
        rs = rs + ' numPhysCells = %s ' % str(self.numPhysCells)

        return rs

class WxArray:
    r"""WxArray(base, wxdata, arrayName, dumpNo, flavor) -> WxArray

    Reads array-data of array with qualified named ``arrayName`` from
    dump ``dumpNo`` for a run called ``base``. The array type can be
    one of 'numeric' for Numeric arrays or 'numpy' for Numpy arrays.
    The variable ``wxdata`` is a reference to the core WxData object.
    This object should only be created using the `WxData.read` method.

    WxArray objects behave like python arrays and hence support the
    python slice syntax. Note that only on array access is any data
    read from the HDF5 file."""

    def __init__(self, wxdata, base, arrayName, frame, flavor, comboSolver):
        self.base = base
        self.arrayName = arrayName
        self.frame = frame
        self.wxdata = wxdata
        # get pointer to data from open file handle
        fq = arrayName
        dps = 'self.wxdata.fh.root.%s' % comboSolver
        dps = dps + '.__getattr__("%s")' % fq
        self.dp = eval(dps)
        #self.dp.flavor = flavor
        self.name = self.dp._v_name
        self.fullShape = self.dp.shape
        self.numComponents = self.fullShape[-1]
        self.shape = self.fullShape[:-1]
        self.onGrid = self.dp._v_attrs.vsMesh

        self.grid = wxdata.grid(self.onGrid, comboSolver)

    def load(self):
        r"""load() -> array

        Loads the complete array data into memory. Should not be used
        if the array is too big to fit into memory.
        """
        return self.dp.read()

    def __getitem__(self, key):
        r"""Slice array and read data slice into memory
        """
        
        # ensure that the dimensions of the slice are correct
        if len(key) != len(self.fullShape):
            raise TypeError

        # return actual data
        return self.dp.__getitem__(key)

    def __repr__(self):
        r"""Array represtation for printing.
        """
        return self.dp.__repr__()

def Pn(n, x):
    p0 = 1.0
    p1 = x

    if n==0: return p0
    if n==1: return p1
    # initiliaze recurrence
    pn = 0.0
    pn1 = p1
    pn2 = p0
    for i in range(2,n+1):
        # use recurrence relation to compute P_n
        pn = (x*(2.*i-1.)*pn1 - (i-1.)*pn2)/(1.*i)
        pn2 = pn1
        pn1 = pn
    return pn

def compQuadPoints(n):
    x = []
    d = 2.0/n
    x.append(-1.0+0.5*d)
    for i in range(1,n):
        x.append(x[i-1]+d)
    return x

class WxDGArray:
    r"""WxDGArray(WxArray : arr, int : meqn) -> WxDGArray

    Returns an object representing an interpolated version of a
    WxArray to quadrature points. An optional grid object can be
    provided so the solution is interpolated on that grid.
    """

    def __init__(self, arr, meqn, grid=None):
        self.onGrid = arr.onGrid # name of grid on which array lives
        self.grid = arr.grid # grid object
        
        ndims = arr.grid.ndims
        # compute spatial order
        if ndims == 1:
            so = self.spatialOrder = arr.numComponents/(meqn*ndims)
        elif ndims == 2:
            so = self.spatialOrder = int(math.sqrt(arr.numComponents/meqn))
        elif ndims == 3:
            so = self.spatialOrder = int(round(pow(arr.numComponents/meqn,1./3.)))

        # construct Legendre polynomials at quadrature points
        xo = compQuadPoints(so)
        legpol = numpy.zeros((so, so), numpy.float)
        for m in range(so):
            for cc in range(so):
                legpol[m,cc] = Pn(m, xo[cc])

        # compute X coordinate of each quadrature point
        Xc = arr.grid.X # cell center X-coordinate
        dx = arr.grid.dx[0] # cell spacing in X
        nx = Xc.shape[0]
        self.Xq = numpy.zeros((so*nx,), numpy.float)

        for j in range(so):
            for i in range(nx):
                self.Xq[so*i+j] = Xc[i] + 0.5*dx*xo[j]

        self.grid.X = self.Xq # reset the X axis
        self.shape = [nx*so]

        # compute Y coordinate of each quadrature point
        if ndims==2:
            Yc = arr.grid.Y # cell center Y-coordinate
            dy = arr.grid.dx[1] # cell spacing in Y
            ny = Yc.shape[0]
            self.Yq = numpy.zeros((so*ny,), numpy.float)

            for j in range(so):
                for i in range(ny):
                    self.Yq[so*i+j] = Yc[i] + 0.5*dy*xo[j]

            self.grid.Y = self.Yq # reset the Y axis
            self.shape.append(ny*so)

        # compute Z coordinate of each quadrature point
        if ndims==3:
            Yc = arr.grid.Y # cell center Y-coordinate
            dy = arr.grid.dx[1] # cell spacing in Y
            ny = Yc.shape[0]
            self.Yq = numpy.zeros((so*ny,), numpy.float)
            
            for j in range(so):
                for i in range(ny):
                    self.Yq[so*i+j] = Yc[i] + 0.5*dy*xo[j]
                    
            self.grid.Y = self.Yq # reset the Y axis
            self.shape.append(ny*so)
            
            Zc = arr.grid.Z # cell center Y-coordinate
            dz = arr.grid.dx[2] # cell spacing in Y
            nz = Zc.shape[0]
            self.Zq = numpy.zeros((so*nz,), numpy.float)

            for j in range(so):
                for i in range(nz):
                    self.Zq[so*i+j] = Zc[i] + 0.5*dz*xo[j]

            self.grid.Z = self.Zq # reset the Y axis
            self.shape.append(nz*so)                

        self.shape.append(meqn)

        nxso = nx*so
        if ndims==1:
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso,meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cc in range(so):
                    qt = numpy.zeros((nx,), numpy.float)
                    for m in range(so):
                        qt = qt + legpol[m][cc]*arr[:,idx.index2(me,m)]
                    qres[cc:nxso:so,me] = qt

        if ndims==2:
            nyso = ny*so
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso, nyso, meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cx in range(so):
                    for cy in range(so):
                        qt = numpy.zeros((nx,ny), numpy.float)
                        for m in range(so):
                            for n in range(so):
                                qt = qt + legpol[m][cx]*legpol[n][cy]*arr[:,:,idx.index3(me,m,n)]
                        qres[cx:nxso:so,cy:nyso:so,me] = qt

        if ndims==3:
            nyso = ny*so
            nzso = nz*so
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso, nyso, nzso, meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so, 0,so, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cx in range(so):
                    for cy in range(so):
                        for cz in range(so):
                            qt = numpy.zeros((nx,ny,nz), numpy.float)
                            for m in range(so):
                                for n in range(so):
                                    for p in range(so):
                                        qt = qt + legpol[m][cx]*legpol[n][cy]*legpol[p][cz]*arr[:,:,:,idx.index4(me,m,n,p)]
                                        
                            qres[cx:nxso:so,cy:nyso:so,cz:nzso:so,me] = qt

                        
        self.res = qres        

    def load(self):
        r"""load() -> array

        Loads the complete array data into memory. Should not be used
        if the array is too big to fit into memory.
        """
        return self.res

    def __getitem__(self, key):
        r"""Slice array and read data slice into memory
        """
        # return data
        return self.res.__getitem__(key)

class WxData:
    r"""WxData(base : string, frm : int, flavor : string) -> WxData

    Provides an interface to read data from a WarpX hyperbolic solver
    simulation with base name ``base`` and frame ``frm``.  Optionally
    the array ``flavor`` can be specified to select the kind of array
    to use (one of numpy or numeric).
    """
    
    def __init__(self, base, frm, flavor='numpy'):
        self.base = base
        self.frame = frm
        self.flavor = flavor
        fn = base + "_%d.h5" % frm

        # ensure file exist
        if not os.path.exists(fn):
            raise "WxData::__init__ : Dump %d of run %s not exist" % (frm, base)

        self.fh = tables.openFile(fn, "r")
        # read in simulation time
        self.time = float(self.fh.root.timeData._v_attrs.time)

    def close(self):
        r"""close() -> None

        Closes the file
        """
        self.fh.close()

    def comboSolver(self):
        r"""comboSolver() -> string

        Returns name of the top comboSolver
        """
        for k in self.fh.root._v_groups:
            if 'type' in self.fh.root._v_groups[k]._v_attrs:
                if self.fh.root._v_groups[k]._v_attrs.type == 'ComboSolver':
                    return k

    def variables(self, comboSolver = None):
        r"""variables(comboSolver : string) -> string []

        Returns list of variables in simulation. If ``comboSolver`` is
        specified it should be the name of the top comboSolver in the
        simulation.
        """

        # get hold of comboSolver
        if comboSolver:
            cs = comboSolver
        else:
            cs = self.comboSolver()

        grp = self.fh.root
        try:
            grp = grp._v_groups[cs]
        except:
            raise "WxData::variables : Group %s does not exist" % cs

        # read all groups which are variables
        vrbls = []
        for c in grp._v_children:
            child = grp._v_children[c]
            if 'vsType' in child._v_attrs:
                if child._v_attrs.__getattribute__('vsType') == 'variable':
                    vrbls.append(c)
        return vrbls
        
    def grids(self, comboSolver = None):
        r"""grids(comboSolver : string) -> string []

        Returns list of grids in simulation.  If ``comboSolver`` is
        specified it should be the name of the top comboSolver in the
        simulation.
        """

        # get hold of comboSolver
        if comboSolver:
            cs = comboSolver
        else:
            cs = self.comboSolver()

        grp = self.fh.root
        try:
            grp = grp._v_groups[cs]
        except:
            raise "WxData::grids : Group %s does not exist" % cs

        ds = []
        for c in grp._v_children:
            child = grp._v_children[c]
            if 'vsType' in child._v_attrs:
                if child._v_attrs.__getattribute__('vsType') == 'mesh':
                    ds.append(c)
        return ds        

    def grid(self, gridName, comboSolver = None):
        r"""grid(gridName : string) -> WxGrid

        Returns a grid object given the name ``gridName``.  If
        ``comboSolver`` is specified it should be the name of the top
        comboSolver in the simulation.
        """

        # get hold of comboSolver
        if comboSolver:
            cs = comboSolver
        else:
            cs = self.comboSolver()

        grp = self.fh.root
        try:
            grp = grp._v_groups[cs]
        except:
            raise "WxData::grids : Group %s does not exist" % cs

        for c in grp._v_children:
            if c == gridName:
                child = grp._v_children[c]                
                if child._v_attrs.__getattribute__('vsType') == 'mesh':
                    break
                else:
                    raise "Object %s exists but is not a grid" % gridName
                    
        lowerBounds = child._v_attrs.vsLowerBounds
        upperBounds = child._v_attrs.vsUpperBounds
        numPhysCells = child._v_attrs.vsNumCells

        return WxGrid(lowerBounds, upperBounds, numPhysCells)

    def read(self, name, comboSolver = None):
        r"""read(name : string) -> WxArray

        Read an array from the output file.  If ``comboSolver`` is
        specified it should be the name of the top comboSolver in the
        simulation.
        """

        if comboSolver:
            return WxArray(self, self.base, name, self.frame, self.flavor, comboSolver)
        else:
            return WxArray(self, self.base, name, self.frame, self.flavor, self.comboSolver())

    def readDG(self, name, meqn, comboSolver = None):
        r"""readDG(name : string, meqn : int) -> WxDGArray

        Read an array from the output file for a discontinuous
        galerkin simulation with ``meqn`` number of equations.  If
        ``comboSolver`` is specified it should be the name of the top
        comboSolver in the simulation.
        """

        wxa = self.read(name, comboSolver)
        return WxDGArray(wxa, meqn)

    def readDGOnGrid(self, name, meqn, grid, comboSolver = None):
        r""" readDGOnGrid(name : string, meqn : int, grid) -> WxDGArray

        Read an array from output file for a discontinous galerkin
        simulations with ``meqn`` number of equations. The grid on
        which the solution should be interpolated should be specified
        in ``grid``.
        """

        wxa = self.read(name, comboSolver)
        return WxDGArray2(wxa, meqn, grid)

    def readDGQuad(self, name, meqn, quad, comboSolver = None):
        r""" readDGOnGrid(name : string, meqn : int, grid) -> WxDGArray

        Read an array from output file for a discontinous galerkin
        simulations with ``meqn`` number of equations. The grid on
        which the solution should be interpolated should be specified
        in ``grid``.
        """

        wxa = self.read(name, comboSolver)
        return WxDGArrayQuad(wxa, meqn, quad)


class WxDGArrayQuad:
    r"""WxDGArray(WxArray : arr, int : meqn, int : quad) -> WxDGArray

    Returns an object representing an interpolated version of a
    WxArray to quadrature points. The number of quadrature points is
    specified so that the solution can be project on that many points.
    """

    def __init__(self, arr, meqn, quad, grid=None):
        self.onGrid = arr.onGrid # name of grid on which array lives
        self.grid = arr.grid # grid object
        
        ndims = arr.grid.ndims
        # compute spatial order
        if ndims == 1:
            so = self.spatialOrder = arr.numComponents/(meqn*ndims)
        elif ndims == 2:
            so = self.spatialOrder = int(math.sqrt(arr.numComponents/meqn))
        elif ndims == 3:
            so = self.spatialOrder = int(round(pow(arr.numComponents/meqn,1./3.)))

        # construct Legendre polynomials at quadrature points
        xo = compQuadPoints(quad)
        legpol = numpy.zeros((so,quad), numpy.float)
        for m in range(so):
            for cc in range(quad):
                legpol[m,cc] = Pn(m, xo[cc])

        # compute X coordinate of each quadrature point
        Xc = arr.grid.X # cell center X-coordinate
        dx = arr.grid.dx[0] # cell spacing in X
        nx = Xc.shape[0]
        self.Xq = numpy.zeros((quad*nx,), numpy.float)

        for j in range(quad):
            for i in range(nx):
                self.Xq[quad*i+j] = Xc[i] + 0.5*dx*xo[j]

        self.grid.X = self.Xq # reset the X axis
        self.shape = [nx*quad]

        # compute Y coordinate of each quadrature point
        if ndims==2:
            Yc = arr.grid.Y # cell center Y-coordinate
            dy = arr.grid.dx[1] # cell spacing in Y
            ny = Yc.shape[0]
            self.Yq = numpy.zeros((quad*ny,), numpy.float)

            for j in range(quad):
                for i in range(ny):
                    self.Yq[quad*i+j] = Yc[i] + 0.5*dy*xo[j]

            self.grid.Y = self.Yq # reset the Y axis
            self.shape.append(ny*quad)

        # compute Y coordinate of each quadrature point
        if ndims==3:
            Yc = arr.grid.Y # cell center Y-coordinate
            dy = arr.grid.dx[1] # cell spacing in Y
            ny = Yc.shape[0]
            self.Yq = numpy.zeros((quad*ny,), numpy.float)

            for j in range(quad):
                for i in range(ny):
                    self.Yq[quad*i+j] = Yc[i] + 0.5*dy*xo[j]

            self.grid.Y = self.Yq # reset the Y axis
            self.shape.append(ny*quad)
            
            Zc = arr.grid.Z # cell center Y-coordinate
            dz = arr.grid.dx[2] # cell spacing in Y
            nz = Zc.shape[0]
            self.Zq = numpy.zeros((quad*nz,), numpy.float)

            for j in range(quad):
                for i in range(nz):
                    self.Zq[quad*i+j] = Zc[i] + 0.5*dz*xo[j]

            self.grid.Z = self.Zq # reset the Y axis
            self.shape.append(nz*quad)                

        self.shape.append(meqn)

        nxso = nx*quad
        if ndims==1:
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso,meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cc in range(quad):
                    qt = numpy.zeros((nx,), numpy.float)
                    for m in range(so):
                        qt = qt + legpol[m][cc]*arr[:,idx.index2(me,m)]
                    qres[cc:nxso:quad,me] = qt

        if ndims==2:
            nyso = ny*quad
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso, nyso, meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cx in range(quad):
                    for cy in range(quad):
                        qt = numpy.zeros((nx,ny), numpy.float)
                        for m in range(so):
                            for n in range(so):
                                qt = qt + legpol[m][cx]*legpol[n][cy]*arr[:,:,idx.index3(me,m,n)]
                        qres[cx:nxso:quad,cy:nyso:quad,me] = qt

        if ndims==3:
            nyso = ny*quad
            nzso = nz*quad
            # allocate memory for interpolated data
            qres = numpy.zeros((nxso, nyso, nzso, meqn), numpy.float)
            # interpolate solution to each quadrature point
            rng = WxRange( (0,meqn, 0,so, 0,so, 0,so) ) # range object
            idx = WxColIndexer(rng)
            for me in range(meqn):
                for cx in range(quad):
                    for cy in range(quad):
                        for cz in range(quad):
                            qt = numpy.zeros((nx,ny,nz), numpy.float)
                            for m in range(so):
                                for n in range(so):
                                    for p in range(so):
                                        qt = qt + legpol[m][cx]*legpol[n][cy]*legpol[p][cz]*arr[:,:,:,idx.index4(me,m,n,p)]
                                        
                            qres[cx:nxso:quad,cy:nyso:quad,cz:nzso:quad,me] = qt

                        
        self.res = qres        

    def load(self):
        r"""load() -> array

        Loads the complete array data into memory. Should not be used
        if the array is too big to fit into memory.
        """
        return self.res

    def __getitem__(self, key):
        r"""Slice array and read data slice into memory
        """
        # return data
        return self.res.__getitem__(key)
