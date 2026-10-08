#!/usr/bin/env bash
# Item 05 -- launched detached; safe to re-run (skips finished cells).
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "05_perturb" "scripts/s05_perturb.py"
