#!/usr/bin/env bash
# Item 02 -- launched detached; safe to re-run (skips finished cells).
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "02_Lmap" "scripts/s02_Lmap.py"
