"""The two configurations s36 did not reach, run on their own."""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s36)

REST = {"a=b /5     a=0.2 c=d=0.05 e=0.02": (0.2, 0.2, 0.05, 0.05, 0.02),
        "all equal  a=b=c=d=0.3 e=0.6": (0.3, 0.3, 0.3, 0.3, 0.6)}

print(f"  {'configuration':38s} {'lam_D':>7} {'lam_v':>7} {'lam_I':>7} "
      f"{'p_pred':>7} {'p_fit':>7} {'rel err':>8} {'r2':>7}")
rows = []
for name, shape in REST.items():
    p_fit, r2, s_lo, s_hi = s36.power_law_p(shape)
    if not np.isfinite(p_fit):
        print(f"  {name:38s}  no crossing at this L range"); continue
    lD, lV, lI = s36.rates(shape, s_lo, s_hi)
    p_pred = lI / (2 * lD - lV)
    rel = (p_pred - p_fit) / p_fit
    rows.append(dict(config=name, lam_Delta=lD, lam_v2=lV, lam_Iss=lI,
                     p_pred=p_pred, p_fit=p_fit, rel_err=rel, r2=r2,
                     s_lo=s_lo, s_hi=s_hi))
    print(f"  {name:38s} {lD:7.4f} {lV:7.4f} {lI:7.4f} {p_pred:7.4f} "
          f"{p_fit:7.4f} {rel:+8.3f} {r2:7.4f}", flush=True)
out = HERE.parent / "results" / "factors" / "factors_rest.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\n  wrote {out}")
