

import math
pi = math.pi
XUP = 3.0
TEND = pi
U = 1.0

# -*- shell-script -*-
<warpx>

  ##
  # Basic parameters
  ##
  Run_name = advect # run name
  Simulation = advection # simulation 
  Real = double # precision type: float/double
  # one of debug, info, warning, error, critical, disabled
  Verbosity = debug

  <advection>
    Type = WxSolver
    Kind = comboSolver

    Time = [0.0, TEND] # start, end times
    Out = 40 # no of output files to write
    Dt = 0.04 # initial dt to try
    UseFixedDt = 0
    
    ##
    # Grid defintion
    ##
    <grid>
    Type = WxGridBox
    Dimensions = 1
    Lower = [-pi] # lower coordinates of block
    Upper = [pi] # upper coordinates of block
    Cells = [nx] # no of cells
    PeriodicDirs = [0] # X is periodic
    </grid>

    ##
    # Array definitions
    ##
    <qold>
      Type = WxVariable
      Kind = parArray
      OnGrid = grid
      NumComponents = 1 # number of components
      GhostCells = [2, 2] # no of ghost cells
    </qold>

    <qnew>
      Type = WxVariable
      Kind = parArray
      OnGrid = grid
      NumComponents = 1 # number of components
      GhostCells = [2, 2] # no of ghost cells
    </qnew>

    ##
    # Subsolver definitions
    ## 

    # hyperbolic solver
    <hyperbolicSolver>
      Type = WxSubSolver
      Kind = hyperSubSolver
      Scheme = wavePropagation
      ReadVars = [qold]
      WriteVars = [qnew]
      OnGrid = grid

      Equations = [advectionint] # equation to solve

      # array initializer
      Initialize = [qold, qnew]
      <InitialCondition>
       Type = WxFunction
       Kind = exprFunc # name of initial-condition function

       a = value
       
       exprList = ["pulse"]
                   
       progn = ["pulse = a*sin(x)"]
       
      </InitialCondition>

      # Wave propagation scheme
      <wavePropagation>
        Kind = wave1d # type of scheme

        Cfl = 1.0 # CFL number to use
        Cflm = 1.0 # maximum CFL allowed

        spatialOrder = 2 # spatial order 1: Gudonov, 2: Lax-Wendroff
        sourceSplitting = 0 # 0: No source, 1: Gudonov splitting, 2: Strang splitting
        # one of minmod, superbee, van-leer, monotonized-centered or beam-warming
        limiter = van-leer
      </wavePropagation>
      
      ##
      # Advection equation
      ##
      <advectionint>
       Type = WxHyperbolicEqn
       Kind = advectionEqn # kind of equations
       ux = U # x-component of advection velocity
      </advectionint>
      
    </hyperbolicSolver>

    ##
    # Subsolver steps
    ##

    <lowerBC>
      Type = WxSubSolver
      Kind = gridFixedBc

      OnGrid = grid
      WriteVars = [qold, qnew]
     
      value = 0.0

      direction = 0
      edge = lower
    </lowerBC>

    <upperBC>
      Type = WxSubSolver
      Kind = bcCopy
      OnGrid = grid
      direction = 0
      edge = upper
    </upperBC>

    # solve hyperbolic equations

    # Set qnew back to qold so it can be used in hyperbolicSolve
    <copier>
      Type = WxSubSolver
      Kind = linearCombiner
      OnGrid = grid
      ReadVars = [qnew]
      coeffs = [1.0]
      WriteVars = [qold]
    </copier>

    ##
    # Subsolver steps
    ##
    # define boundary condition step
    <applyBC>
      Type = WxSubSolverStep

      SubSolvers = [lowerBC, upperBC] # subsolver to apply
    </applyBC>
    

    <solveHyperEqn>
      Type = WxSubSolverStep
      DtFrac = 1.0
      SubSolvers = [hyperbolicSolver]
    </solveHyperEqn>

    <copy>
      Type = WxSubSolverStep
      SubSolvers = [copier]
    </copy>
    
    ##
    # Solver sequence
    ##
    <SolverSequence>
      Type = WxSolverSequence
      PerStep = [solveHyperEqn, copy] # apply at each step
    </SolverSequence>

  </advection>

</warpx>
