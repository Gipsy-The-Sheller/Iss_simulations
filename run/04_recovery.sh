#!/usr/bin/env bash
# Recovery half of Item 3 only -- the alignments are already on disk, so this
# re-measures P_true without re-simulating.  ML is the expensive part.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "04_recovery" "scripts/s04_zone.py" --recovery-only
