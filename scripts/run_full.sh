#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
export R_LIBS_USER="$project_root/.Rlibrary${R_LIBS_USER:+:$R_LIBS_USER}"
exec .venv/bin/python run.py "$@"
