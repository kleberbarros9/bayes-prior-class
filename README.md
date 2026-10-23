# Full reproduction of the Bayesian prior-generating class article

This repository runs the **final analyses**, in English, for the three results sections of *A Bayesian Prior-Generating Class: Construction, Induced Posteriors, and Applications*.

There is **no dataset screening stage and no metrics-only/raw execution mode**. Default execution performs the full experiments, posterior fits, predictive calculations, tables, figures, and the independent precision check. The fixed applications are the eight Roever cohorts at 80 mg/m² and the NCCTG lung cohort. P-spline-PH is the flexible comparator in the final nine-model survival comparison.

## 1. Obtain the project

Clone the authors' repository:

```bash
git clone https://github.com/kleberbarros9/bayes-prior-class.git
cd bayes-prior-class
```

If this release was supplied as a ZIP before the files were committed to GitHub, extract it first and use its `bayes-prior-class` folder. The repository was inspected at commit `7bee27ffb9579e5cc689d9b1c30d440b6436ded9`, when it contained only a README. The supplied release adds the implementation; the GitHub remote is not modified automatically.

Keep the directory structure intact. Commands below are run from the project root. Scripts also resolve their paths correctly when invoked from another working directory.

## 2. Install Linux system dependencies

Use a 64-bit Linux computer with Python 3.10–3.13, R 4.3 or newer, a C++17 compiler, and GNU make. The documented container uses **R 4.5.0**. A practical starting point is 16 GB RAM, with more memory useful for Stan compilation and independent integration. Reserve several GB for posterior checkpoints and predictive outputs; full archives can be much larger than compact summaries.

On Debian/Ubuntu, install the system components once:

```bash
sudo apt-get update
sudo apt-get install -y \
  r-base r-base-dev build-essential gfortran \
  python3 python3-venv python3-pip \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libtbb-dev \
  libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev zlib1g-dev
```

`sudo` is needed only for system packages. Do not run the analysis as root. On a managed cluster, ask the administrator for these components or load the corresponding modules.

Check the basic executables:

```bash
Rscript --version
python3 --version
g++ --version
make --version
```

## 3. Create the project environments

With internet access, run:

```bash
bash scripts/setup_linux.sh
```

This creates `.venv`, installs the pinned Python dependencies, installs the specified R package versions in `.Rlibrary`, checks the C++ toolchain with an isolated toy Stan model, and checks R/Stan syntax, data rules and numerical identities. RStan compiles its models locally; CmdStan and a separate Stan download are not required.

The R installer restores `renv.lock` if present. On the first installation, it resolves the specified root versions and writes a lockfile for the complete installed dependency set. **Commit that generated lockfile** when publishing the verified environment. This release does not pretend that its root-version list is already a complete transitive lockfile. `installed_R_packages.csv` and per-analysis `sessionInfo.txt` record the actual environment.

Root R versions: RcppParallel 5.1.7, rstan 2.32.7, loo 2.8.0, posterior 1.6.1, ggplot2 3.5.1, gridExtra 2.3, splines2 0.5.4, digest 0.6.37, jsonlite 1.9.1, survival 3.8-3. See `R/install_dependencies.R`. Local validation can also use compatible preinstalled packages; a different environment must be reported when interpreting numerical replication.

## 4. Run all three full results sections

Foreground execution:

```bash
bash scripts/run_full.sh --workers 4
```

For execution that continues after the SSH terminal disconnects:

```bash
bash scripts/start_full.sh --workers 4
```

Monitor the main log:

```bash
tail -n 40 -F results/full_run.log
```

Or open the checkpoint progress display:

```bash
bash scripts/monitor_progress.sh
```

`Ctrl+C` closes the log viewer or monitor; it does not stop a detached analysis. Progress bars count completed chains/stages and do not assume that all stages have equal cost.

`--workers` limits concurrently launched survival chains (maximum four per model) and the meta-analysis chain-core setting. BLAS/OpenMP threads are restricted to one per process. The three sections run sequentially, keeping the workflow within a predictable resource budget. Computation time depends strongly on hardware and diagnostic retries; validation is not an estimate of the full runtime.

## 5. Run an individual full section

```bash
bash scripts/run_full.sh --section simulation --workers 4
bash scripts/run_full.sh --section meta-analysis --workers 4
bash scripts/run_full.sh --section survival-model --workers 4
```

Every command uses the full scientific configuration. `--validate` is a separate code/data check that does **not** produce or replace the manuscript experiments:

```bash
bash scripts/run_full.sh --validate
```

### Simulation

`R/run_simulation.R` runs 16 Poisson scenarios with 5,000 replications each, two-baseline SBC with 5,000 replications per baseline and 100 posterior draws per replication, 12 independent-censoring scenarios with 1,000 replications, and the special-case rejection-sampling/benchmark experiments. `python/conflict.py` adds the 30 directed upper-tail conflict configurations with 2,000 common replications per scenario. Completed major simulation stages are checkpointed.

The simulation's base seed is 20261019; its Poisson, SBC, censoring and rejection stages use offsets +1, +2, +3 and +4, respectively. The conflict experiment uses 20261024. Shared datasets are used for paired method comparisons. Quantile matching does not imply robustness under upper-tail conflict.

### Meta-analysis

`R/run_meta_analysis.R` reads the frozen full `dat.roever2022` snapshot, fixes the dose at 80 mg/m², and verifies eight distinct study cohorts. It fits the historical MAP models, computes QROG and exponential-odds priors plus six conventional/borrowing comparators, completes leave-one-study-out and patient-level calculations, generates the new-study predictive block, samples 40,000 posterior values for conditional interval checks, and fits the joint exponential-odds pooled model.

Full historical fits start at four chains of 4,000 iterations, with 2,000 warm-up, and allow diagnostic retries. The base seed is 20261020 with deterministic fold offsets. Numerical quadrature/refinement and simulation uncertainty are recorded separately. A failed diagnostic is not silently accepted. Valid completed folds and posterior blocks are reused.

### Survival models

`R/prepare_survival.R` starts from all 228 lung records, excludes the one missing ECOG record, and constructs the 227-patient/164-death analysis from scratch. It centers/scales age and ECOG and encodes female sex. The stored likelihood unit is a quarter; calculations use the corresponding event-density Jacobian, while plots and predictive follow-up use years. The calibration target is `H0(1 year)`, with median 1 and 95th percentile 8, including the induced AFT location priors and their Jacobians.

The nine full models are Quadratic-Weibull, QROG-Weibull, Weibull-PH, Exponential-PH, Lognormal-AFT, Loglogistic-AFT, Piecewise-RW1, Mspline-PH and Bayesian P-spline-PH. Initial full fits use four chains of 16,000 iterations, with 8,000 warm-up per chain. Retries increase iterations; diagnostics require rank Rhat ≤ 1.01, bulk/tail ESS ≥ 4,000, no divergences or depth hits, and BFMI ≥ 0.3. Results remain incomplete if these criteria are not met.

The P-spline baseline is a cubic seven-function log-hazard expansion, normalized at one year, with identified/noncentered RW2 smoothing, proper linear-component and log-smoothing-scale priors, and a constant-hazard tail after three years. Its fit is checked by higher-order quadrature and event-time inversion before final scoring. See `R/pspline_prepare.R`, `R/pspline_math.R` and `models/pspline_ph.stan`.

For every model the workflow generates PSIS-LOO, LOOIC, WAIC, parameterization-specific DIC, targeted direct LOO refits if required, parameter and hazard-ratio summaries, standardized survival curves, Cox–Snell residuals, 40,000 replicated cohorts, future-patient predictions, prior/posterior plots, traces and densities. Cohort predictions use the documented stratified reverse-KM censoring approximation; those checks do not establish independent censoring.

`python/refine_precision.py` then compares the three Weibull priors using common-point RQMC at 2^18 and 2^20 points, 12 scrambles per resolution, and selected direct leave-one-out posterior integrations. It reports pointwise/paired uncertainty and resolution stability. Separate-chain rankings can change under this refinement; no outcome is selected or forced by the code. ELPD and LOOIC are the same criterion up to a factor of −2, WAIC estimates a related target, and DIC depends on parameterization.

## 6. Locate the outputs

| Location | Content |
|---|---|
| `results/simulation/` | Full simulation RDS, scenario CSVs, report, SBC and validation PDF |
| `results/simulation/conflict/` | Directed-conflict summaries, paired comparisons and replicate-level table |
| `results/meta_analysis/oncology/` | Full meta-analysis RDS, interval figures/tables, report and session information |
| `results/full_lung/<model>/` | Base-model chains, fit diagnostics, criteria, summaries, predictions and figures |
| `results/pspline_lung/P-spline-PH/` | Full P-spline fit, numerical audits, criteria, predictions and figures |
| `results/survival-model/` | Combined nine-model CSV/LaTeX tables and paired comparisons |
| `results/lung_refinement_v1/precision/` | Independent common-point and direct LOO precision checks |
| `results/logs/` | Stage-specific execution logs |
| `results/workflow_status.json` | Overall completion/error state |
| `prepared/` | Recreated analysis data and Python integration inputs |
| `reference_outputs/` | Previous compact estimates for comparison, never used as fitted input |

Survival PDF figures include standardized survival, four-panel residual/PPC/ELPD checks, prior/posterior baseline hazard, parameter densities, traces and future-patient predictions. Figures are regenerated after all nine curves are available so overlays contain P-spline-PH.

Fresh Stan compilations, dependency versions and Monte Carlo trajectories can change the last reported digits. Check `reference_comparison.csv`, posterior MCSEs, patient-level uncertainty and independent integration, rather than requiring identical rounded rankings.

## 7. Resume safely after interruption

Run the **same command** in the same checkout:

```bash
bash scripts/start_full.sh --workers 4
```

Completed compatible sections, simulation stages, meta-analysis folds, survival chains, predictive chunks and numerical-integration replicates are reused. An interrupted unfinished chain restarts; Stan sampling within that chain is not resumed mid-trajectory. A postprocessing error does not require refitting accepted chains.

A process lock prevents two workflow runners from writing into the same tree. Chain lock files record host/PID; the runner removes a lock only when its recorded process is no longer alive on the same host. If an orphan chain is still running, let it finish before restarting. Do not delete RDS files to repair a plotting failure.

Source/data fingerprints reject incompatible reuse after implementation changes. Preserve the previous `results`, `prepared` and compiled model files in a separate directory before starting a scientifically changed configuration. Previously generated screening/raw checkpoints are not required or automatically imported into this standalone release.

## 8. Export compact and complete results

After execution completes:

```bash
bash scripts/export_results.sh
sha256sum -c exports/SHA256SUMS
```

This produces `exports/results_compact.tar.gz` without RDS or full replicated-cohort matrices, and `exports/results_with_rds.tar.gz` with posterior checkpoints. Both are written to temporary filenames, checked with `tar`, and renamed only after successful completion. The compact archive is intended for sharing tables, figures, diagnostics and logs. It cannot replace posterior draws for an independent posterior reanalysis.

## 9. Offline L200 execution

Datasets are included, so execution needs no internet after dependencies are installed. A virtual environment or compiled R/Stan library copied from another machine is not automatically portable: Linux distribution, architecture, Python/R versions and compiler libraries must be compatible.

On L200 with the required R packages already installed, keep its existing `R_LIBS_USER` or set that variable to the correct cluster library. Use an offline Python wheelhouse built on a **compatible Linux/Python system**:

```bash
# Connected compatible host:
python3 -m pip download -r requirements.txt -d wheelhouse
# Transfer the project plus wheelhouse/ to L200, then there:
python3 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheelhouse -r requirements.txt
bash scripts/run_full.sh --validate
bash scripts/start_full.sh --workers 4
```

R package installation/compilation must already have been completed on L200 or through the administrator's compatible offline package mechanism. Do not run the network-based setup script on an offline host. Stan models compile on L200 from the included `.stan` files.

From the personal computer, the established jump connection is:

```bash
ssh -J kleber@200.17.137.190:10224 kleber@l200
```

Run a standalone copy in its own folder, such as `~/bayes-prior-class-full`; this release does not overwrite the existing `projeto_walfrido_python/current` outputs.

## 10. Optional container and GitHub publication

On an internet-connected machine with Docker:

```bash
docker build -t bayes-prior-class-full .
mkdir -p results exports
docker run --rm -v "$PWD/results:/analysis/results" -v "$PWD/exports:/analysis/exports" \
  bayes-prior-class-full
```

The container installs the declared R/Python dependencies and performs validation during the build. The full experiments start when the container runs. Host directories preserve outputs after container removal.

To publish this release from a checkout of the authors' repository, review and commit the source, data snapshots, generated lockfile and validation records; `.gitignore` excludes large results and checkpoints:

```bash
git status --short
git add README.md SOURCE_PROVENANCE.md VERSION Dockerfile .dockerignore .gitignore \
  R models python scripts tests data reference_outputs requirements.txt run.py
# Add renv.lock after the verified environment has generated it:
[ ! -f renv.lock ] || git add renv.lock
git commit -m "Add full simulation, meta-analysis and survival reproduction workflows"
git push origin HEAD
```

See `tests/VALIDATION.md` for the checks actually completed for this release. A successful syntax/numerical validation is not a claim that all production Monte Carlo experiments have been rerun.
