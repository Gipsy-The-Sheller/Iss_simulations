#!/usr/bin/env bash
# Recovery targets for experiment 6 (MPBoot + in-process Fitch), over the
# alignments already on disk.  Incremental; MPBoot is one subprocess per
# alignment so this is the slow half.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
launch "06b_recovery" "scripts/s06_condindep.py" --recovery
