#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
mkdir -p exports
# Stage each archive, verify it, then atomically expose the completed filename.
metadata=(reference_outputs)
[[ ! -f installed_R_packages.csv ]] || metadata+=(installed_R_packages.csv)
tar --exclude='*.rds' --exclude='posterior_predictive.csv.gz' --exclude='RUN.lock' -czf exports/results_compact.tar.gz.partial results "${metadata[@]}"
 tar -tzf exports/results_compact.tar.gz.partial >/dev/null
mv exports/results_compact.tar.gz.partial exports/results_compact.tar.gz
tar --exclude='RUN.lock' -czf exports/results_with_rds.tar.gz.partial results prepared "${metadata[@]}"
tar -tzf exports/results_with_rds.tar.gz.partial >/dev/null
mv exports/results_with_rds.tar.gz.partial exports/results_with_rds.tar.gz
sha256sum exports/results_compact.tar.gz exports/results_with_rds.tar.gz > exports/SHA256SUMS
ls -lh exports/results_*.tar.gz
