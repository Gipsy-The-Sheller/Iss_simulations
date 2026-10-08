#!/usr/bin/env bash
# Item 6 -- conditional independence: Pr(Y | pbar, Iss) vs Pr(Y | pbar)
# Runs the crossed simulation, then the model comparison, in one background job.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "06_condindep" -c '
import subprocess, sys
for a in ([], ["--analyze"]):
    r = subprocess.run([sys.executable, "scripts/s06_condindep.py"] + a)
    if r.returncode: sys.exit(r.returncode)
'
