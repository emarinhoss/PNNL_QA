

MI = 1.0
ME = value
Q = 1.0
GAMMA = 1.66666666667
LIGHT_SPEED = 1.0
EPS0 = 1.0

# -*- shell-script -*-
<warpx>

  ##
  # Basic parameters
  ##
  Run_name = ssrecon_wv # run name
  Simulation = ssrecon # simulation 
  Real = double # precision type: float/double
  # one of debug, info, warning, error, critical, disabled
  Verbosity = debug

  <ssrecon>
    Type = WxSolver
    Kind = comboSolver

    Time = [0.0, 400.0] # start, end times
    Out = 40 # no of output files to write
    Dt = 0.04 # initial dt to try
    UseFixedDt = 0
    ##
    # Grid defintion
    ##
    <grid>
    Type = WxGridBox
    Dimensions = 2

    Lower = [-12.8, -6.4] # lower coordinates of block
    Upper = [12.8, 6.4] # upper coordinates of block
    Cells = [512, 256] # no of cells

    PeriodicDirs = [0] # X is periodic
    </grid>

    ##
    # Array definitions
    ##
    <qold>
      Type = WxVariable
      Kind = parArray

      OnGrid = grid
      NumComponents = 18 # number of components
      GhostCells = [2, 2] # no of ghost cells
    </qold>

    <qnew>
      Type = WxVariable
      Kind = parArray

      OnGrid = grid
      NumComponents = 18 # number of components
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

      Equations = [eulerElc, eulerIon, maxwell] # equation to solve
      Sources = [lorentzElc, lorentzIon, ionCurrents, elcCurrents, rhoC] # source terms

      # array initializer
      Initialize = [qold, qnew]
      <InitialCondition>
        Type = WxFunction
        Kind = twoFluidRecon # name of initial-condition function

        qe = -Q # electron charge
        qi = Q  # ion charge
        me = ME # electron mass = 1/25
        mi = MI  # ion mass

        gas_gamma = GAMMA # gas constant

      </InitialCondition>

      # Wave propagation scheme
      <wavePropagation>
        Kind = wave2d # type of scheme

        Cfl = 0.45 # CFL number to use
        Cflm = 0.5 # maximum CFL allowed

        spatialOrder = 2 # spatial order 1: Gudonov, 2: Lax-Wendroff
        sourceSplitting = 1 # 1: Gudonov splitting, 2: Strang splitting
        # one of minmod, superbee, van-leer, monotonized-centered or beam-warming
        limiter = van-leer
      </wavePropagation>

      ##
      # Electron equations
      ##
      <eulerElc>
        Type = WxHyperbolicEqn
        Kind = eulerEqn

        gas_gamma = GAMMA
      </eulerElc>

      ##
      # Ion equations
      ##
      <eulerIon>
        Type = WxHyperbolicEqn
        Kind = eulerEqn

        gas_gamma = GAMMA
      </eulerIon>

      ##
      # Maxwell equations
      ##
      <maxwell>
        Type = WxHyperbolicEqn
        Kind = phMaxwellEqn

        c0 = LIGHT_SPEED
        gamma = 1.0 # error propagation speed for div(B)
        chi = 1.0
      </maxwell>

      ##
      # Lorentz source term for electron
      ##
      <lorentzElc>
        Type = WxHyperbolicSrc
        Kind = lorentzForces

        InpRange = [0,1,2,3, 10,11,12,13,14,15]
        OutRange = [1,2,3,4]

        mass = ME
        charge = -Q
      </lorentzElc>

      ##
      # Lorentz source term for ions
      ##
      <lorentzIon>
        Type = WxHyperbolicSrc
        Kind = lorentzForces

        InpRange = [5,6,7,8, 10,11,12,13,14,15]
        OutRange = [6,7,8,9]

        mass = MI
        charge = Q
      </lorentzIon>

      ##
      # Current sources for electric field from electrons
      ##
      <elcCurrents>
        Type = WxHyperbolicSrc
        Kind = currents

        InpRange = [1,2,3]
        OutRange = [10,11,12]

        mass = ME
        charge = -Q
        epsilon0 = EPS0
      </elcCurrents>

      ##
      # Current sources for electric field from ions
      ##
      <ionCurrents>
        Type = WxHyperbolicSrc
        Kind = currents

        InpRange = [6,7,8]
        OutRange = [10,11,12]

        charge = Q
        mass = MI
        epsilon0 = EPS0
      </ionCurrents>

      <rhoC>
        Type = WxHyperbolicSrc
        Kind = chargeDensity

        InpRange = [0,5] # electron and ion density respectively
        OutRange = [16] # phi equation

        qi = Q
        qe = -Q
        mi = MI
        me = ME
        epsilon0 = EPS0
        chi = 1.0
      </rhoC>

    </hyperbolicSolver>


    # Set qnew back to qold so it can be used in hyperbolicSolve
    <copier>
      Type = WxSubSolver
      Kind = linearCombiner

      OnGrid = grid
      ReadVars = [qnew]
      coeffs = [1.0]
      WriteVars = [qold]
    </copier>

   # use q to see if wpe, wpi, wce, wci are resolved
   <chk_freq>
     Type = WxSubSolver
     Kind = checkFreqTimeStepWv

     charge = Q
     ionmass = MI
     elcmass = ME
     epsilon0 = EPS0
     
     OnGrid = grid

     ReadVars = [qold]

   </chk_freq>


    <lowerBC>
      Type = WxSubSolver
      Kind = bcConductingTwoFluid

      OnGrid = grid
      WriteVars = [qold, qnew]

      direction = 1
      edge = lower
    </lowerBC>

    # define lower boundary condition applicator
    <upperBC>
      Type = WxSubSolver
      Kind = bcConductingTwoFluid

      OnGrid = grid
      WriteVars = [qold, qnew]

      direction = 1
      edge = upper
    </upperBC>

    ##
    # Subsolver steps
    ##

    # solve hyperbolic equations
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

    <chkfreq>
      Type = WxSubSolverStep
      SubSolvers = [chk_freq]
    </chkfreq>

    ##
    # Solver sequence
    ##
    <SolverSequence>
      Type = WxSolverSequence
      PerStep = [chkfreq, applyBC, solveHyperEqn, applyBC, copy] # apply at each step

    </SolverSequence>

  </ssrecon>

</warpx>
