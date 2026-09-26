# PNNL_QA — Uncertainty Quantification for WARPX Plasma Simulations

Scripts used at PNNL (2011) to study how uncertainty in physical input
parameters propagates through fluid plasma simulations run with the
WARPX code (the University of Washington multi-fluid code, not the LBNL
particle-in-cell code of the same name).

**Status:** archived research code. It targets Python 2, a 2011 WARPX build
and specific clusters, and has not been run since. See
[Fixes made in 2026](#fixes-made-in-2026) for bugs corrected after the
original runs.

Three ways of sampling the uncertain parameter are compared:

| Method | Abbreviation | Idea |
| ------ | ------------ | ---- |
| Monte Carlo | **MC** | Random samples, all on one grid resolution. |
| Multilevel Monte Carlo | **MMC / MLMC** | Random samples on a hierarchy of grids. Most samples use the coarse grid; fewer samples use finer grids to correct the coarse-grid estimate. |
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
    ├── euler_dispersion/    Driver for the dispersive Euler case, plus convergence and cost plots
    └── reconnection/        Driver for the reconnection case
```

Every script starts with a comment saying what it reads and writes.

### `advection_testcase/`
- `input2.py`: WARPX template for 1D advection. The initial condition is
  `a*sin(x)`, and the amplitude `a` is the uncertain parameter (`value`),
  drawn from U(0.5, 3).
- `input.py`: a 2D Gaussian advection template. No script uses it.
- `preprocessing.py` (MC), `mmc_runs_setup.py` (MMC) and `pcm_runs_setup.py`
  (PCM) each create one folder per sample, containing a `.pin` file (and
  `weight.w` for PCM).
- `run_mc.sh`, `run_cases.sh` (MMC) and `run_pcm.sh` preprocess and run every
  sample folder.

### `dispersion_testcase/`
The same workflow applied to `input3.py`, which is Euler equations with a
dispersive source term. The uncertain parameter is the initial ion velocity
`u_i`, drawn from U(1e-10, 1e-4).
- `input.py` and `input2.py` are copies of the advection templates. No
  script here uses them.
- `exec.sh` runs the full pipeline: setup for all three methods, then all runs.
- `run_mmc{36,72,150,520}.sh` run MMC batches of different sizes, from
  folders renamed by hand to `advect_mmc_<N>_*`.
- `mc_stat_5000.py` resamples a large MC run in `MC_RUNS/` to show how the
  mean and variance converge as the sample count grows.
- `mc_values.txt` and `mmc_values.txt` record the sampled parameter values.
- `elecacoustic-wv.pin` / `.inp` are a standalone deterministic case.

### `reconnection/`
2D two-fluid magnetic reconnection. The uncertain parameter is the
ion-to-electron mass ratio `MI/ME`, drawn from U(25, 100) with `MI = 1`.
One exception: the current MC script, `cases_setup/preprocessing.py`, fixes
`MI/ME = 25` and varies the speed of light instead.
- `cases_setup/`: templates (`input.py`, `input2.py`), setup scripts for
  MC, MMC and PCM, PBS (`cray*.qsub`) and Moab (`batch_pcm.msub`) job
  scripts, and helpers:
  - `run.sh` and `check.sh` check whether runs finished.
  - `qchange.py` cancels a range of batch jobs.
- `post_process/`:
  - `calc_flux.py` runs inside a single run folder. It integrates |B_y| along
    the mid-plane to get the reconnected flux for each output frame, and
    writes `<prefix>.dat`. With no arguments it processes every level or PCM
    run it finds in the folder.
  - `stats_*.py` and `flux-vis-*.py` gather those files across run folders
    and compute the mean and variance. Each takes an optional glob for the
    run folders (see [Run folder names](#run-folder-names)).
  - `mass_ratio_plots.py` compares the three methods.
- `early_study/`: the first version of this study.
  - `preprocessing.m`, `postprocessing.m`, `restart_postprocessing.m` and
    `compare_plot.m` are MATLAB scripts that vary both the mass ratio and the
    speed of light.
  - `meanflux_and_variation.py` is a Python version of the post-processing.
  - The folder also holds the job scripts used at that stage.

### `mlmc_driver/`
Adaptive MLMC in Python (`multiMonteCarlo.py` and `original_mmc_warpx.py`),
following Giles (2008), *Multilevel Monte Carlo path simulation*, Operations
Research 56(3). Level `l` uses `nx * M**l` cells, and `P_l` is the quantity
of interest on that grid. Starting from level 0, the driver:
1. Runs `N` samples on the new level. Each sample also runs level `l-1`
   with the same parameter, so it can form the correction
   `Y_l = P_l - P_{l-1}`.
2. Estimates the variance `V_l` of `Y_l` on each level and sets the number of
   samples per level to

   `N_l = ceil( 2 / e**2 * sqrt(V_l / C_l) * sum_k sqrt(V_k * C_k) )`

   where `C_l = (nx * M**l)**gamma` is the cost of one sample on level `l`.
   `gamma` is 2 for the 1D Euler case and 3 for the 2D reconnection case.
3. Adds samples where needed and tests convergence with
   `max(|mean(Y_{L-1})| / M, |mean(Y_L)|) < (M - 1) * e / sqrt(2)`.
   If that test fails, it adds level `L+1` and goes back to step 1.

The MLMC estimate of the mean is `sum_l mean(Y_l)`. In the reconnection
driver, each iteration writes the running sums (`suml1.dat` … `suml5.dat`),
the sampled values (`rand_vals.dat`) and the settings (`run_info.dat`).

`euler_dispersion/Nl_and_cost_plot.py` and
`euler_dispersion/log_varn_and_log_mean.py` produce the usual MLMC
diagnostics:
- samples per level
- MLMC cost against standard MC cost, as a function of tolerance
- decay of the mean and variance with level

`reconnection/batchjob.msub` runs the reconnection driver on a Moab cluster.

## Run folder names

Setup scripts put the sampled value in each folder name. A numeric prefix
identifies the method:

| Case | MC | MMC | PCM |
| ---- | -- | --- | --- |
| advection, dispersion | `advect_002_U_<value>` | `advect_001_U_<value>` | `advect_003_U_<value>` |
| reconnection (setup scripts) | `recon_008_MR_<MI/ME>_c0_<c>` | `recon_006_MR_<MI/ME>_c0_<c>` | `recon_001_ME_<ME>` |
| reconnection (post-processing defaults) | `recon_001*` | `recon_001*` | `recon_002*` |

The reconnection runs were renamed by hand between setup and
post-processing, so the prefixes don't match. Pass the right glob to the
post-processing scripts, for example `python stats_mc.py "./recon_008_MR_*"`.

MMC grid sizes (levels 0 / 1 / 2):

| Case | Level 0 | Level 1 | Level 2 |
| ---- | ------- | ------- | ------- |
| advection | 100 | 200 | 400 |
| dispersion | 100 | 500 | 1000 |
| reconnection | 400×200 | 800×400 | 1600×800 |

All samples run on level 0, the first half also on level 1, and the first
quarter also on level 2.

## Output files

| File | Written by | Contents |
| ---- | ---------- | -------- |
| `*.pin` | setup scripts | Sampled parameter and grid size followed by the WARPX template; `$wxpp` turns it into `*.inp`. |
| `weight.w` | PCM setup | Quadrature weight of that collocation node. The weights sum to 1. |
| `<prefix>_<frame>.h5` | WARPX | Solution at each output frame; read with `wxdata.py`. |
| `<prefix>.dat` | `calc_flux.py` | Reconnected flux per frame, normalised to 0.2 at t = 0. |
| `mc_mean.dat`, `mc_vari.dat` | `stats_mc.py` | MC mean and unbiased sample variance per frame. |
| `pcm_mean_2.dat`, `pcm_vari_2.dat` | `stats_pcm.py` | PCM mean and variance per frame. |
| `mean.dat`, `vari.dat` | `stats_mmc.py` | MLMC mean; three variance estimates as columns. |
| `mc_stat_<N>.dat` | `mc_stat_5000.py` | Mean and variance per grid cell from N resampled MC runs. |

## Requirements

- **WARPX** executables, referenced through environment variables:
  - `$wxpp`: preprocessor that turns a `.pin` file into a `.inp` file.
  - `$warpx`, `$warpxs`, `$warpxser`, `$warpxp`: parallel, serial and MPI
    builds of the solver.
- **Python 2** with `numpy` (older than 1.24, which removed `numpy.float`),
  `matplotlib` (`pylab`) and PyTables (`tables`). The scripts use Python 2
  syntax, such as `print x`.
- `wxdata.py`, the WARPX HDF5 reader. It is in `mlmc_driver/*/` and must be
  on your `PYTHONPATH` for the post-processing scripts.
- MATLAB, only for `reconnection/early_study/*.m`. Those scripts need the
  `load_data_new` reader, which is not in this repo.
- A PBS or Moab scheduler, for the `*.qsub` and `*.msub` job scripts. These
  contain the original author's account and email; edit them before
  submitting.

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

The MC and MMC setup scripts draw random samples. To make a run
reproducible, set `SEED` at the top of the script to an integer.

For reconnection:
1. Run the setup and job scripts in `reconnection/cases_setup/`.
2. In each finished run folder, run `python calc_flux.py`.
3. From the directory that holds the run folders, run the `stats_*.py` and
   `flux-vis-*.py` scripts.

Generated run folders and WARPX outputs (`*.h5`, `*.log`) are ignored by git.

## Fixes made in 2026

These changes correct bugs in the original scripts. Results produced in 2011
with the old versions may differ.

- **MLMC sample counts** (`mlmc_driver/*/multiMonteCarlo.py`,
  `Nl_and_cost_plot.py`): the optimal `N_l` used `sum sqrt(V_l / n_l)`
  instead of Giles's `sum sqrt(V_l * C_l)`, so it took far fewer samples
  than the tolerance requires. It now also uses a per-sample cost model
  (`gamma`).
- **MLMC variance output** (`mlmc_driver/reconnection/multiMonteCarlo.py`):
  `varn.dat` contained the mean. The variance also used `1/(N**2-N)` with
  integer `N`, which is 0 in Python 2.
- **MLMC cost plot** (`Nl_and_cost_plot.py`):
  - The MC cost now uses the variance of `P_L`, not of the corrections.
  - The y axis now shows `e**2 * cost`, as labelled.
- **MLMC diagnostic scripts** (`Nl_and_cost_plot.py`,
  `log_varn_and_log_mean.py`): these called `preprocess()` with an argument
  list that no version of `original_mmc_warpx.py` accepts. They were written
  for the Euler driver, so they now live in `euler_dispersion/`.
  `log_varn_and_log_mean.py` also had the same integer-division bug.
- **`stats_mmc.py`**: one variance loop indexed with the wrong variable, so
  it counted level 2 twice and level 1 not at all.
- **`flux-vis-pcm.py`**: the error bars showed `sqrt(E[x^2])` instead of the
  standard deviation.
- **`mc_stat_5000.py`**:
  - It crashed because `os` was never imported.
  - It never reset its sums between sample sizes.
  - It divided by the wrong count.
  - It never saved its results.
- **Early-study error bars** (`*.m`, `meanflux_and_variation.py`): these
  plotted the variance instead of the standard deviation.
- **MMC setup scripts**: levels 1 and 2 now get exactly half and a quarter
  of the samples. Previously they got one extra each.
- **Setup scripts**: output files were never explicitly closed (`out.close`
  without `()`).
- **Cleanup**:
  - `replace_text.py` no longer fails with an indentation error.
  - The dead draft `Ommc_warpx.py` was removed.
  - Six copies of the flux calculation were merged into `calc_flux.py`.
