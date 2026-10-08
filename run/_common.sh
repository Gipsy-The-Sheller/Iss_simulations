#!/usr/bin/env bash
# Shared preamble for the nohup launchers.
#   usage:  bash run/03_iss_pdist.sh
# Logs:  logs/<name>.log     PID:  logs/<name>.pid
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO" || exit 1

PY="${PY:-python}"
export PYTHONIOENCODING=utf-8
export PYTHONUNBUFFERED=1
# All simulation goes through AliSim inside IQ-TREE; point at it explicitly.
export IQTREE="${IQTREE:-E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/iqtree3/iqtree3.exe}"

mkdir -p logs results

launch() {                     # launch <name> <scripts/...py>
  local name="$1"; shift
  local log="logs/${name}.log" pid="logs/${name}.pid"
  if [[ -f "$pid" ]] && kill -0 "$(cat "$pid")" 2>/dev/null; then
    echo "[$name] already running (pid $(cat "$pid")) -- not starting again"
    return 1
  fi
  nohup "$PY" "$@" >>"$log" 2>&1 &
  echo $! >"$pid"
  echo "[$name] started pid $(cat "$pid")  -> $log"
}
