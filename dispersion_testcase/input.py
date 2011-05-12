


XUP = 20.0
YUP = 3.0
V = 0.0 # Y-advection speed
TEND = 10.0

<warpx>
##
# Basic parameters
##
 Run_name = rname # run name
 Simulation = comboSolver # simulation 
 Real = double # precision type: float/double
 # one of debug, info, warning, error, critical, disabled
 Verbosity = debug

 <comboSolver>
   Type = WxSolver
   Kind = comboSolver   
    
   Time = [0.0, TEND]
   Dt = 0.1*TEND
   Out = 40
 
   # grid on which to solve equations
   <grid>
     Type = WxGridBox

     Lower = [0.0, 0.0]
     Upper = [XUP, YUP]
     Cells = [nx, ny]

     PeriodicDirs = [0,1]
   </grid>

   # arrays for storing solution
   <qnew>
     Type = WxVariable
     Kind = parArray

     OnGrid = grid
     NumComponents = 1
     GhostCells = [2, 2]
   </qnew>
   
   <q>
     Type = WxVariable
     Kind = parArray

     OnGrid = grid
     NumComponents = 1
     GhostCells = [2, 2]
   </q>

   # define the hyperbolic subsolver
   <hyperbolic>
     Type = WxSubSolver
     Kind = hyperSubSolver

     OnGrid = grid
     ReadVars = [q]
     WriteVars = [qnew]

     Scheme = wavePropagation
     Equations = [advection]

     # define function for initial conditions
     Initialize = [q, qnew] # arrays to initialize
     <InitialCondition>
       Type = WxFunction
       Kind = exprFunc # name of initial-condition function

       # constants for use in expressions
       xc = XUP/2.0
       yc = YUP/2.0

       # list of preiliminary expressions to execute
       progn = ["r2 = (x-xc)^2 + (y-yc)^2"]
       # list of expressions: one per component to initiliaize
       exprList = ["exp(-10*r2)"]

     </InitialCondition>

     # define scheme parameters
     <wavePropagation>
       Type = WxHyperbolicScheme
       Kind = wave2d # type of scheme

       Cfl = 0.49 # CFL number to use
       Cflm = 0.5 # maximum CFL allowed

       spatialOrder = 2 # spatial order 1: Gudonov, 2: Lax-Wendroff
       sourceSplitting = 0 # 0: no source, 1: Gudonov splitting, 2: Strang splitting
       # one of no-limiter minmod, superbee, van-leer, monotonized-centered or beam-warming
       limiter = no-limiter
     </wavePropagation>     

     # define equation parameters
     <advection>
       Type = WxHyperbolicEqn
       Kind = advectionEqn # kind of equations

       ux = U # x-velocity
       uy = V # y-velocity
       
     </advection>

   </hyperbolic>

   # define subsolver to copy qold to q
   <copier>
     Type = WxSubSolver
     Kind = linearCombiner

     OnGrid = grid

     ReadVars = [qnew]
     coeffs   = [1.0]
     WriteVars = [q]

   </copier>

   # define homogenous step
   <homogeneous>
     Type = WxSubSolverStep

     DtFrac = 1.0 # time fraction to apply this update step
     SubSolvers = [hyperbolic] # subsolver to apply
   </homogeneous>

   # define copy step
   <copy>
     Type = WxSubSolverStep

     SubSolvers = [copier] # subsolver to apply
   </copy>

   # define subsolver sequence to apply
   <SolverSequence>
     Type = WxSolverSequence

     PerStep = [homogeneous, copy] # sequence of steps
   </SolverSequence>

 </comboSolver>

</warpx>
