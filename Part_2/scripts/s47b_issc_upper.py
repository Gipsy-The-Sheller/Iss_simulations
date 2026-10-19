"""Iss.c recomputed on the *falling* crossing, which is Xia's definition.

Two things were wrong in s47.

First, the crossing.  s47 took the first crossing of 95% -- the tree length at
which recovery climbs up through it.  That is not what the index means.  Xia
et al. (2003) define it as the point at which the sequences *begin to fail*:

    "the critical Iss value (referred to hereafter as Iss:c) at which the
     sequences will begin to fail to recover the true tree"

Recovery against tree length rises and then falls: at large tree lengths all
three topologies approach the same saturated distribution and the signal that
distinguishes them vanishes.  So the level is crossed twice, and Iss.c is the
*second* crossing, the largest tree length still recoverable.  Reading the
first one instead inverts the dependence on SeqLen -- s47 gave Iss.c falling
with SeqLen, contradicting both Xia's Fig. 3 and the theory, in which

    Delta(s) sqrt(L) = z sqrt(v2(s))

has two roots precisely because Delta(s) is not monotone.

Second, the crossing was read off 50 replicates pointwise.  Recovery
plateaus, and the plateau sits within about two standard errors of the level,
so the crossing landed wherever the noise happened to fall and the resulting
Iss.c was not even monotone in L.  The falling limb is therefore fitted --
a logistic in log tree length, which is the shape the fall has -- and the
crossing is read off the fit rather than off a single point.

Xia's fallback for the asymmetrical topology is implemented too, since on this
tree recovery at 500 bp never reaches 95% at all:

    "For the asymmetrical topology ... there is often no tree length at which
     the true tree is recovered 100%.  However, there is always a tree length
     at which Ptrue is the largest.  When Ptrue decreases to 95% of the
     maximum Ptrue value, the tree length is taken as TLC."

Both Xia's switching rule and the uniform 95%-of-max rule are reported.  The
switching rule is not one measure, because which branch of it applies changes
with L, so the uniform column is the one to read.

Outputs
    results/fragment_pool/issc_upper.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
OUT = PART2 / "results" / "fragment_pool"


def crossing(sub, level):
    """Iss where recovery falls through `level`, and the scale there.

    Pointwise, from the last grid scale still at or above the level.  A
    fitted falling limb was tried first and is not usable: the plateau sits
    within about two standard errors of the level, so the fit is free to put
    the crossing almost anywhere and it extrapolated to scale 1e8.  The
    pointwise answer is noise-limited but at least lands on the grid.

    The precision of this number is the replicate count, not the method.  s47
    used 50; the endpoint values are stable and the intermediate ones are not.
    """
    s = sub.sort_values("scale")
    iss, rec = s.iss.values, s.recovery.values
    idx = np.flatnonzero(rec >= level)
    if len(idx) == 0:
        return None, None
    i = idx[-1]
    if i + 1 >= len(iss) or rec[i] == rec[i + 1]:
        return None, None
    w = (level - rec[i]) / (rec[i + 1] - rec[i])
    return (float(iss[i] + w * (iss[i + 1] - iss[i])),
            float(s.scale.values[i] + w * (s.scale.values[i + 1]
                                           - s.scale.values[i])))


def main() -> None:
    d = pd.read_csv(OUT / "issc_curve.csv")
    print("Iss.c on the fitted falling crossing (Xia's 'begin to fail')\n")
    print(f"  {'L':>7} {'max rec':>8} {'level':>7} {'Iss.c':>8} {'scale':>8}")
    rows = []
    for L in sorted(d.L.unique()):
        sub = d[d.L == L]
        pmax = float(sub.recovery.max())
        lam = 0.95 if pmax >= 0.95 else 0.95 * pmax
        c_uni, sc_uni = crossing(sub, 0.95 * pmax)
        c_xia, sc_xia = crossing(sub, lam)
        rows.append(dict(L=L, max_rec=pmax, level=lam,
                         iss_c_uniform=c_uni, scale_uniform=sc_uni,
                         iss_c_xia=c_xia, scale_xia=sc_xia))
        fmt = lambda v: "--" if v is None else f"{v:.4f}"            # noqa: E731
        print(f"  {L:7d} {pmax:8.3f} {lam:7.3f} {fmt(c_uni):>8} "
              f"{'--' if sc_uni is None else f'{sc_uni:8.3f}'}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "issc_upper.csv", index=False)

    got = out.dropna(subset=["iss_c_uniform"]).sort_values("L")
    print("\n  Iss.c(L), uniform 95%-of-max rule:")
    for r in got.itertuples():
        print(f"    L={int(r.L):>7}   Iss.c = {r.iss_c_uniform:.4f}"
              f"   at scale {r.scale_uniform:.3f}"
              f"  (tree length {2.564 * r.scale_uniform:.2f})")
    if len(got) >= 2:
        inc = bool((got.iss_c_uniform.diff().dropna() > 0).all())
        print(f"    -> {'INCREASING' if inc else 'NOT monotone'} in L")
        if not inc:
            print("       Xia's Fig. 3 has it increasing; a non-monotone")
            print("       result here means the fit is still noise-limited,")
            print("       not that the index behaves differently.")

    a = out[out.L == 500].iss_c_uniform.iloc[0]
    b = out[out.L == 10_000].iss_c_uniform.iloc[0]
    if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
        meta = pd.read_csv(OUT / "meta.csv")
        if b > a:
            m = (meta.iss >= a) & (meta.iss < b)
            print(f"\n  the two calibrations disagree on Iss in "
                  f"({a:.4f}, {b:.4f}):")
            print(f"  {m.mean():.4f} of the pool, {int(m.sum())} fragments, "
                  f"{int(m.sum()) * 500} bp")
        else:
            print(f"\n  Iss.c(500)={a:.4f} >= Iss.c(10k)={b:.4f}: the "
                  "calibrations are ordered")
            print("  the other way on this tree, so there is no band where "
                  "they disagree.")


if __name__ == "__main__":
    main()
