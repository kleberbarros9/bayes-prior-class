#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
mkdir -p results
nohup bash scripts/run_full.sh "$@" >> results/full_run.log 2>&1 < /dev/null &
echo "$!" > results/launcher.pid
printf 'Launcher PID: %s\nLog: tail -n 40 -F %s/results/full_run.log\n' "$!" "$project_root"
