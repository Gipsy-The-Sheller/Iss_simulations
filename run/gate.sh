#!/usr/bin/env bash
# G0 gate -- every downstream run assumes this passes.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "g0_gate" "scripts/s00_gate.py"
