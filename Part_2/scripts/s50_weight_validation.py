"""Is the vote table in s49 what CASTER-site actually does?

s49 predicts CASTER's recovery probability on a quartet in closed form, from
an assumed weight table.  That prediction is falsifiable: run the binary on
simulated quartets and compare.  If the table is wrong the curves separate
immediately, and the w3 result of s49 would be built on a fiction.

The test is run on the Farris-zone quartet, where the margin at w3 = -1/2 is
only +0.0012 per site.  A wrong table would almost certainly put that margin
on the other side of zero, and the binary would then get the topology *wrong*
at long alignments instead of right -- a sign, not a shade.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "caster_weights"
WORK = PART2 / "_caster" / "weights"
CASTER = Path("E:/Codes/ASTER/exe/caster-site.exe")

_spec = importlib.util.spec_from_file_location("s49", str(HERE / "s49_caster_weights.py"))
s49 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s49)
m = s49.m

_s45 = importlib.util.spec_from_file_location("s45", str(HERE / "s45_caster_discard.py"))
s45 = importlib.util.module_from_spec(_s45)
_s45.loader.exec_module(s45)

BASES = "ACGT"
TRUE_SPLIT = "{A,B}|{C,D}"
L_GRID = (10_000, 30_000, 100_000, 300_000, 1_000_000)
N_REP = 60


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    p0 = m.pattern_batch(s49.FARRIS, m.SPLITS["12|34"])[0]
    V = s49.vote_table(-0.5)
    rng = np.random.default_rng(20261016)

    print("Does the closed form in s49 reproduce the CASTER-site binary?")
    print(f"  Farris quartet {s49.FARRIS}, JC, w3 = -1/2, {N_REP} replicates\n")
    print(f"  {'L':>9} {'predicted':>10} {'observed':>9} {'se':>7} {'diff':>8}")

    rows = []
    for L in L_GRID:
        pred = s49.power(V, p0, L)
        ok = 0
        for r in range(N_REP):
            idx = rng.choice(256, size=L, p=p0)
            d = np.stack([(idx >> 6) & 3, (idx >> 4) & 3, (idx >> 2) & 3, idx & 3])
            fa = WORK / f"v_{L}_{r}.fa"
            with open(fa, "w") as f:
                for t, n in enumerate("ABCD"):
                    f.write(f">{n}\n{''.join(BASES[c] for c in d[t])}\n")
            lst = WORK / f"v_{L}_{r}_list.txt"
            lst.write_text(str(fa.resolve()) + "\n")
            tre = WORK / f"v_{L}_{r}.nw"
            if tre.exists():
                tre.unlink()
            subprocess.run([str(CASTER), "-i", str(lst), "-o", str(tre)],
                           capture_output=True, text=True, cwd=str(WORK))
            got = None
            if tre.exists() and tre.read_text().strip():
                got = s45.recovered_split(tre.read_text().strip().splitlines()[-1])
            ok += int(got == TRUE_SPLIT)
        obs = ok / N_REP
        se = np.sqrt(obs * (1 - obs) / N_REP)
        rows.append(dict(L=L, predicted=pred, observed=obs, n=N_REP, se=float(se)))
        print(f"  {L:9d} {pred:10.4f} {obs:9.4f} {se:7.4f} "
              f"{obs - pred:+8.4f}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "validation.csv", index=False)

    d = pd.DataFrame(rows)
    z = (d.observed - d.predicted) / d.se.clip(lower=1e-9)
    print(f"\n  mean |z| = {np.abs(z).mean():.2f} over {len(d)} points")
    print("  The prediction is exact and has no fitted parameter, so a")
    print("  systematic offset would mean the vote table is not CASTER's.")


if __name__ == "__main__":
    main()
