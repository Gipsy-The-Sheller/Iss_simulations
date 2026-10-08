#!/usr/bin/env bash
# Launch everything in order, sequentially inside one background job.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "all" -c '
import subprocess, sys
for s in ["scripts/s00_gate.py","scripts/s01_design.py","scripts/s02_Lmap.py",
          "scripts/s03_iss_pdist.py","scripts/s04_zone.py","scripts/s05_perturb.py"]:
    print("=== "+s, flush=True)
    subprocess.run([sys.executable, s], check=False)
'
