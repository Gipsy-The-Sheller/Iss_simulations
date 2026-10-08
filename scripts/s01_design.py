"""Item 0 -- uniform-Iss sampling design (Part II Theorem 11).

No simulation here.  Computes mu(s) exactly, the density p_s(s) = mu'(s)/Z,
the inverse-CDF sampler, and the marginal-flatness diagnostic (exact mu vs
second-order mu).  Writes results/design/density.csv and design.json, and the
s-grids that s02-s04 consume.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import design as D                                       # noqa: E402
from issbench.models import build_Q                                    # noqa: E402

OUT = ROOT / "results" / "design"
OUT.mkdir(parents=True, exist_ok=True)

CONFIGS = [
    dict(name="q4_sym_homog",   N=4, fi=0.35, alpha=np.inf, s_min=0.05, s_max=60.0),
    dict(name="q4_sym_a1",      N=4, fi=0.35, alpha=1.0,  s_min=0.05, s_max=60.0),
    dict(name="q4_short_homog", N=4, fi=0.10, alpha=np.inf, s_min=0.05, s_max=40.0),
    dict(name="q4_short_a05",   N=4, fi=0.10, alpha=0.5,  s_min=0.05, s_max=200.0),
]


def main():
    Q = build_Q("JC69")
    pi = np.ones(4) / 4
    meta = {}
    rows = []
    for cfg in CONFIGS:
        d = D.quartet_design(Q, pi, cfg["fi"], N=cfg["N"], alpha=cfg["alpha"],
                             s_min=cfg["s_min"], s_max=cfg["s_max"])
        ss, mus = d.mu_grid(240)
        aff = np.array([D.mu_affine(d, s) for s in ss])
        dens = d.density(ss)
        for s, m, a, p in zip(ss, mus, aff, dens):
            rows.append(dict(config=cfg["name"], s=s, mu_exact=m,
                             mu_second_order=a, p_s=p))
        # target grid: 60 points uniform in Iss between 10% and 90% of the range
        lo, hi = mus.min(), mus.max()
        targets = lo + np.linspace(0.10, 0.90, 60) * (hi - lo)
        ss_sorted = np.argsort(mus)
        grid = np.interp(targets, mus[ss_sorted], ss[ss_sorted])
        meta[cfg["name"]] = dict(
            fi=cfg["fi"], alpha=("inf" if np.isinf(cfg["alpha"]) else cfg["alpha"]),
            s_min=cfg["s_min"], s_max=cfg["s_max"],
            mu_min=float(lo), mu_max=float(hi), mu_span=float(hi - lo),
            max_second_order_error=float(np.abs(mus - aff).max()),
            s_grid=[float(x) for x in grid],
            iss_grid=[float(x) for x in targets],
        )
        print(f"{cfg['name']:16s} mu in [{lo:.4f},{hi:.4f}]  "
              f"max|2nd-order err|={np.abs(mus-aff).max():.4f}  "
              f"gap to exact ceiling {1-0.6620-hi:+.4f}", flush=True)

    pd.DataFrame(rows).to_csv(OUT / "density.csv", index=False)
    (OUT / "design.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT/'density.csv'} and design.json")


if __name__ == "__main__":
    main()
