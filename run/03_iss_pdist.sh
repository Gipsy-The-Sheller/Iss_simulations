#!/usr/bin/env bash
# Item 03 -- launched detached; safe to re-run (skips finished cells).
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "03_iss_pdist" "scripts/s03_iss_pdist.py"
