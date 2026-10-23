#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
command -v Rscript >/dev/null || { echo 'Install R first; see README.md, Linux system dependencies.' >&2; exit 1; }
command -v python3 >/dev/null || { echo 'Python 3 is required.' >&2; exit 1; }
python3 -c 'import sys; assert (3,10) <= sys.version_info[:2] <= (3,13), "Use Python 3.10–3.13 for the pinned dependency set."'
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
export R_LIBS_USER="$project_root/.Rlibrary"
Rscript R/install_dependencies.R "$project_root"
Rscript R/check_toolchain.R
.venv/bin/python run.py --validate
printf '\nSetup complete. Start: bash scripts/run_full.sh\n'
