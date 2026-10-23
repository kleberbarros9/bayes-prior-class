#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
if [[ "${1:-}" == "--once" ]]; then exec .venv/bin/python python/monitor.py; fi
watch -n 10 "$project_root/.venv/bin/python $project_root/python/monitor.py"
