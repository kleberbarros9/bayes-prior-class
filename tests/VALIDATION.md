# Validation completed for this release

The following checks were executed while preparing the package:

| Check | Outcome |
|---|---|
| Python source compilation and CLI argument parsing | Passed |
| Bash syntax for every shell script | Passed |
| R parsing for all R files | Passed |
| Stan syntax translation for both survival model files | Passed |
| Actual C++ compilation and isolated Stan sampling installation check | Passed; the tiny check is not a posterior diagnostic assessment |
| Public input rules | Passed: eight fixed-dose cohorts; 227 complete lung patients and 164 deaths |
| Generated-prior normalization | Passed in R and Python |
| Common-target AFT Jacobians against numerical derivatives | Passed |
| P-spline quadrature, zero-spline hazard and inverse identities, including zero and tail times | Passed |
| Simulation workflow on its internal test configuration | Passed; test outputs are separate from full experiments |
| Repeating that simulation workflow | Passed; major-stage checkpoints reused with identical inferential summaries |
| Directed conflict experiment at the full specified replication count | Completed: 30 method/scenario configurations, 2,000 replications per scenario |
| Conflict coverage and RMSE against the supplied manuscript outputs | Maximum absolute discrepancy 0 for both quantities |
| Independent integration API on a small validation grid | Passed; does not certify production precision |
| Top-level `run.py --validate` entry point | Passed |

The local validation runtime used R 4.3.3 and compatible preinstalled packages; `validation_R_package_versions.csv` identifies them. The recommended production installer specifies the archived analysis root-package versions and the container uses R 4.5.0. No claim is made that the Docker build or a new clean installation of every pinned dependency was executed here.

The complete production simulation, all historical meta-analysis chains and all nine full survival fits were **not rerun** during packaging. Logs distinguish numerical/code checks from production inference. The reference output files come from the supplied completed runs and are never substituted as new analysis inputs.

The simulation checkpoint check invokes the source's internal reduced test configuration solely for software validation. The distributed workflow has no screening, raw, or reduced scientific execution profile: its section commands always run the full configuration.

The RcppParallel/TBB layout of this temporary local runtime required a locally built RcppParallel 5.1.7 for its compiler check. The project installer explicitly pins this package; temporary runtime relocations are not included in the project.
