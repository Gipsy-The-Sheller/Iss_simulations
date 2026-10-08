#!/usr/bin/env bash
# Item 04 -- launched detached; safe to re-run (skips finished cells).
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "04_zone" "scripts/s04_zone.py"
