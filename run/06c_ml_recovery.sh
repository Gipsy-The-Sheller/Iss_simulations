#!/usr/bin/env bash
# Item 6, third pass -- ML recovery (IQ-TREE) over the E6 alignments.
# Sharded across SHARDS processes, then merged.  Incremental per shard, so
# re-running costs only what is missing.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

SHARDS="${SHARDS:-14}"
export SHARDS PY
LOG="logs/06c_ml_recovery.log"
PID="logs/06c_ml_recovery.pid"

if [[ -f "$PID" ]] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "[06c_ml_recovery] already running (pid $(cat "$PID")) -- not starting again"
  exit 1
fi

nohup bash -c '
  pids=""
  for i in $(seq 0 $((SHARDS-1))); do
    "$PY" scripts/s06_ml_recovery.py --shard "$i" --of "$SHARDS" &
    pids="$pids $!"
  done
  for p in $pids; do wait "$p"; done
  "$PY" scripts/s06_ml_recovery.py --merge --of "$SHARDS"
' >>"$LOG" 2>&1 &

echo $! >"$PID"
echo "[06c_ml_recovery] started pid $(cat $PID) with $SHARDS shards -> $LOG"
