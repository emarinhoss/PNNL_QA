# PNNL_QA — Uncertainty Quantification for WARPX Plasma Simulations

Scripts used at PNNL (2011) to study how uncertainty in physical input
parameters propagates through fluid plasma simulations run with the
WARPX code (the University of Washington multi-fluid code, not the LBNL
particle-in-cell code of the same name).

Three ways of sampling the uncertain parameter are compared:

| Method | Abbreviation | Idea |
| ------ | ------------ | ---- |
| Monte Carlo | **MC** | Random samples, all on one grid resolution. |
| Multilevel Monte Carlo | **MMC / MLMC** | Random samples on a hierarchy of grids (e.g. 100 / 200 / 400 cells). Most samples use the coarse grid; fewer samples use finer grids to correct the coarse-grid estimate. |
| Probabilistic Collocation | **PCM** | Deterministic Clenshaw–Curtis nodes (`pointsN`) and weights (`weightsN`), with N = 33 or 65. |

Each method produces the mean and variance of a quantity of interest, such as
the solution profile or the reconnected magnetic flux.

## Repository layout

```
.
├── advection_testcase/      1D linear advection, random amplitude (MC / MMC / PCM)
├── dispersion_testcase/     1D Euler with a dispersive (electron-acoustic) source term
├── reconnection/            Magnetic reconnection (GEM-type) with uncertain mass ratio
│   ├── cases_setup/         Generate run folders and submit jobs (MC / MMC / PCM)
│   ├── post_process/        Compute reconnected flux, statistics and plots
│   └── early_study/         First study (MATLAB + older job scripts): varies mass ratio and speed of light
└── mlmc_driver/             Adaptive MLMC driver: runs WARPX itself and chooses levels and sample counts
    ├── euler_dispersion/    Driver for the dispersive Euler case
    └── reconnection/        Driver for the reconnection case, plus convergence and cost plots
```

### `advection_testcase/`
- `input2.py`: WARPX template for 1D advection. The initial condition is
  `a*sin(x)`, and the amplitude `a` is the uncertain parameter (`value`),
  drawn from U(0.5, 3).
- `input.py`: a 2D Gaussian advection template. The setup scripts do not
  use it.
- `preprocessing.py` (MC, `advect_002_*`), `mmc_runs_setup.py` (MMC,
  `advect_001_*`), `pcm_runs_setup.py` (PCM, `advect_003_*`): each creates
  one folder per sample, containing a `.pin` file (and `weight.w` for PCM).
- `run_mc.sh`, `run_cases.sh` (MMC), `run_pcm.sh`: preprocess and run every
  sample folder.

### `dispersion_testcase/`
The same workflow applied to `input3.py`, which is Euler equations with a
dispersive source term. The uncertain parameter is the initial ion velocity
`u_i`, drawn from U(1e-10, 1e-4).
- `exec.sh` runs the full pipeline: setup for all three methods, then all runs.
- `run_mmc{36,72,150,520}.sh` run MMC batches of different sizes.
- `mc_stat_5000.py` resamples a large MC run to study how the error
  converges as the sample count grows.
- `mc_values.txt` and `mmc_values.txt` record the sampled parameter values.

### `reconnection/`
2D two-fluid magnetic reconnection. The uncertain parameter is the
ion-to-electron mass ratio `MI/ME`, drawn from U(25, 100) with `MI = 1`.
One exception: the current MC script, `cases_setup/preprocessing.py`, fixes
`MI/ME = 25` and varies the speed of light instead.
- `cases_setup/`: templates (`input.py`, `input2.py`), setup scripts for
  MC, MMC and PCM, PBS (`cray*.qsub`) and Moab (`batch_pcm.msub`) job
  scripts, and helpers that check whether runs finished (`run.sh`,
  `check.sh`).
- `post_process/`:
  - `calc_flux_*.py` run inside a single run folder. They integrate |B_y| on
    the mid-plane to get the reconnected flux, and write it to a `.dat`
    file.
  - `stats_*.py` and `flux-vis-*.py` gather those files across run folders
    and compute the mean and variance.
  - `mass_ratio_plots.py` compares the methods.
- `early_study/`: the first version of this study. `preprocessing.m`,
  `postprocessing.m` and `restart_postprocessing.m` are MATLAB scripts that
  also vary the speed of light. The folder also holds
  `meanflux_and_variation.py` and the job scripts used at that stage.

### `mlmc_driver/`
Adaptive MLMC in Python (`multiMonteCarlo.py` and `original_mmc_warpx.py`).
Starting from level 0, the driver:
1. Runs N samples on each level, with a grid refinement factor of M.
2. Estimates the variance per level and computes the optimal number of
   samples N_l.
3. Adds samples where needed and tests convergence. If the estimate has not
   converged, it adds a finer level.

`reconnection/Nl_and_cost_plot.py` and `reconnection/log_varn_and_log_mean.py`
produce the usual MLMC diagnostics: samples per level, cost against
tolerance, and mean and variance decay with level.
`reconnection/batchjob.msub` runs the driver on a Moab cluster.

## Requirements

- **WARPX** executables, referenced through environment variables:
  - `$wxpp`: preprocessor that turns a `.pin` file into a `.inp` file.
  - `$warpx`, `$warpxs`, `$warpxser`, `$warpxp`: parallel, serial and MPI
    builds of the solver.
- **Python 2** with `numpy`, `matplotlib` (`pylab`) and PyTables (`tables`).
  The scripts use Python 2 syntax, such as `print x`.
- `wxdata.py`, the WARPX HDF5 reader. It is in `mlmc_driver/*/` and must be
  on your `PYTHONPATH` for the post-processing scripts.
- MATLAB, only for `reconnection/early_study/*.m`. Those scripts need the
  `load_data_new` reader, which is not in this repo.
- A PBS or Moab scheduler, for the `*.qsub` and `*.msub` job scripts.

## Typical workflow

Run all scripts from inside their own directory. They read templates and
quadrature files from the current working directory.

```bash
cd advection_testcase
python preprocessing.py      # MC:  creates advect_002_U_* folders
python mmc_runs_setup.py     # MMC: creates advect_001_U_* folders (levels 0/1/2)
python pcm_runs_setup.py     # PCM: creates advect_003_U_* folders + weight.w
sh run_mc.sh; sh run_cases.sh; sh run_pcm.sh
```

For reconnection, copy the scripts you need from `reconnection/post_process/`
into the directory that holds the `recon_*` run folders, then run them there.

Generated run folders and WARPX outputs (`*.h5`, `*.log`) are ignored by git.
