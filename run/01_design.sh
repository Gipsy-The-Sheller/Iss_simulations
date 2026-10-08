#!/usr/bin/env bash
# Item 01 -- launched detached; safe to re-run (skips finished cells).
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "01_design" "scripts/s01_design.py"
