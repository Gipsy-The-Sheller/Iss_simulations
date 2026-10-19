"""Figures for the completed experiments.

House style mirrors ref2.png: serif panel titles, sans-serif axis labels, muted
palette, thin grey spines, light grid, dashed reference lines, cream shading for
the calibration region.

    python scripts/s07_figures.py            # all figures whose data is ready
    python scripts/s07_figures.py fig_zone   # one
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.ticker import LogLocator                             # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench.results import res                                  # noqa: E402

RES = ROOT / "results"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

# ------------------------------------------------------------------ style
C = dict(black="#1a1a1a", red="#B22222", orange="#E8871E", blue="#3B7EA1",
         green="#2E8B57", purple="#7B4F9D", grey="#8C8C8C")
CREAM = "#FDF3E0"


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        # mathtext defaults to DejaVu Sans whatever font.sans-serif says, which
        # would leave the maths inside an Arial label in a different face.  Make
        # the maths Arial too; Arial has Greek, but no blackboard bold, so the
        # one \mathbb in the figures falls back rather than erroring.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        "mathtext.tt": "Arial",
        "mathtext.cal": "Arial:italic",
        "mathtext.default": "it",
        "mathtext.fallback": "stix",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#4d4d4d",
        "axes.labelcolor": "#1a1a1a",
        "axes.grid": True,
        "grid.color": "#dcdcdc",
        "grid.linewidth": 0.6,
        "xtick.color": "#4d4d4d",
        "ytick.color": "#4d4d4d",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#cccccc",
        "figure.dpi": 130,
        "savefig.dpi": 200,
    })


def panel(ax, s):
    """Panel letter, top-left inside the axes.

    Descriptive headings belong in the caption, not on the figure: systematic
    biology figures carry no title of their own, and a panel heading would have
    to be repeated in the caption anyway.
    """
    # above the axes, not inside: the top-left corner is where the legends go
    ax.set_title(s, loc="left", fontsize=11, fontweight="bold",
                 family="serif", color="#1a1a1a", pad=5)


def title(ax, s):
    """Draw only the panel letter of a descriptive heading.

    Call sites still pass the old headings ("A. Location vs L"); the prose has
    moved into the figure captions, so only the label before the first full
    stop is drawn -- not the first character, which would turn the two-digit
    grid labels A1/A2/B1/B2 into A/A/B/B.
    """
    panel(ax, s.strip().split(".", 1)[0])


def save(fig, name):
    """Fixed canvas: an out-of-bounds artist must not be allowed to inflate the
    saved image (bbox_inches='tight' did exactly that here)."""
    p = FIGS / f"{name}.png"
    try:
        fig.tight_layout()
    except Exception:                                                 # noqa: BLE001
        pass
    fig.savefig(p, bbox_inches=None)
    plt.close(fig)
    print(f"  wrote {p.relative_to(ROOT)}")


# ================================================================ Fig 0
def fig_design():
    f = res(RES / "design", "density.csv")
    if not f.exists():
        return
    d = pd.read_csv(f)
    order = ["q4_sym_homog", "q4_sym_a1", "q4_short_homog", "q4_short_a05"]
    labels = {"q4_sym_homog": "quartet, $f_i$=0.35, homogeneous",
              "q4_sym_a1": "quartet, $f_i$=0.35, $\\Gamma$ $\\alpha$=1",
              "q4_short_homog": "quartet, $f_i$=0.10, homogeneous",
              "q4_short_a05": "quartet, $f_i$=0.10, $\\Gamma$ $\\alpha$=0.5"}
    cols = [C["black"], C["red"], C["blue"], C["orange"]]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7))
    for c, lab, col in zip(order, [labels[k] for k in order], cols):
        s = d[d.config == c].sort_values("s")
        axes[0].plot(s.s, s.mu_exact, color=col, lw=1.8, label=lab)
        axes[0].plot(s.s, s.mu_second_order, color=col, lw=1.0, ls=":")
        axes[1].plot(s.s, s.p_s, color=col, lw=1.8)
        # density INDUCED on I_ss: p_s(s)/|dmu/ds|, which is flat iff the
        # sampler really returns a uniform I_ss.  Plotting p_s against mu (as
        # this panel used to) does not test that at all.
        # both curves are masked where dmu/ds degenerates: mu_exact is stored
        # rounded, so beyond saturation the numerical gradient is exactly 0 and
        # the quotient is infinite
        lab = ("exact", "2nd order") if c == order[0] else (None, None)
        for mu, ls, lw, lb in ((s.mu_exact, "-", 1.8, lab[0]),
                               (s.mu_second_order, ":", 1.0, lab[1])):
            g = np.gradient(mu.to_numpy(), s.s.to_numpy())
            den = s.p_s.to_numpy() / np.abs(g)
            ok = (np.abs(g) > 0.02) & (den < 3.0)
            axes[2].plot(mu[ok], den[ok], color=col, lw=lw, ls=ls, label=lb)
    axes[0].set_xscale("log")
    axes[0].axhline(0.661995, color=C["grey"], ls="--", lw=1)
    axes[0].set_xlabel("scale $s$"); axes[0].set_ylabel(r"$\mu(s)=\mathbb{E}[I_{ss}]$")
    title(axes[0], "A. Design curve $\\mu(s)$ (solid exact, dotted 2nd order)")
    axes[0].legend(fontsize=6.5, loc="upper left")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("scale $s$"); axes[1].set_ylabel(r"$p_s(s)=\mu'(s)/Z$")
    title(axes[1], "B. Sampling density")
    axes[2].set_xlabel(r"$\mu(s)$")
    axes[2].set_ylabel(r"density on $I_{ss}$")
    axes[2].set_ylim(0, 2.4)
    title(axes[2], "C. Induced density on $I_{ss}$ (flat = uniform)")
    axes[2].legend(fontsize=6.5, loc="upper left")
    save(fig, "fig0_design_density")


# ================================================================ Fig 1
def fig_lmap():
    f = res(RES / "exp1_Lmap", "cells.csv")
    if not f.exists():
        return
    d = pd.read_csv(f)
    d = d[d.n_seq == 4]
    stats = {"median": ("$I_{ss}$ (median)", "YlOrBr_r"),
             "std": ("SD$(I_{ss})$", "PuBu_r")}
    # one colourbar per ROW on a dedicated narrow column: attaching a bar to
    # each of the four axes put the right column's y-label underneath the left
    # column's colourbar
    fig = plt.figure(figsize=(9.8, 6.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.04],
                          wspace=0.12, hspace=0.44, right=0.93)
    for r, (skey, (slab, cmap)) in enumerate(stats.items()):
        im = None
        for c, shape in enumerate(("sym", "asym")):
            sub = d[d["cell"].str.startswith(shape)]
            sub = sub.assign(TL=sub["cell"].str.extract(r"TL([\d.]+)").astype(float))
            piv = sub.pivot_table(index="L", columns="TL", values="iss_A",
                                  aggfunc=skey)
            ax = fig.add_subplot(gs[r, c])
            im = ax.imshow(piv.values, aspect="auto", origin="lower", cmap=cmap)
            ax.set_xticks(range(len(piv.columns)))
            ax.set_xticklabels([f"{x:g}" for x in piv.columns],
                               rotation=90, fontsize=6)
            ax.set_yticks(range(len(piv.index)))
            if c:
                ax.set_yticklabels([])         # shares the row's y-scale
            else:
                ax.set_yticklabels([f"{int(x):d}" for x in piv.index], fontsize=7)
                ax.set_ylabel("sequence length")
            ax.grid(False)
            title(ax, f"{'AB'[r]}{c+1}. {slab}, "
                      f"{'equal' if shape=='sym' else 'graded'} pendants")
            ax.set_xlabel("total tree length")
        fig.colorbar(im, cax=fig.add_subplot(gs[r, 2]))
    save(fig, "fig1_Lmap_heatmaps")

    # ---- profiles: read the invariances off directly
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))
    sub = d[d["cell"].str.startswith("sym")]
    sub = sub.assign(TL=sub["cell"].str.extract(r"TL([\d.]+)").astype(float))
    tls = sorted(sub.TL.unique())
    tls = [tls[3], tls[7], tls[11], tls[16], tls[22]]   # from the grid
    cols = [C["black"], C["red"], C["orange"], C["blue"], C["green"]]
    for tl, col in zip(tls, cols):
        s = sub[np.isclose(sub.TL, tl)]
        med = s.groupby("L").iss_A.median()          # location: median, not the
        q1 = s.groupby("L").iss_A.quantile(0.25)     # 200 raw replicates drawn
        q3 = s.groupby("L").iss_A.quantile(0.75)     # as a line
        axes[0].plot(med.index, med.values, "o-", ms=3.5, color=col,
                     label=f"TL={tl:g}")
        axes[0].fill_between(med.index, q1.values, q3.values, color=col,
                             alpha=0.18, lw=0)
        sd = s.groupby("L").iss_A.std()
        axes[1].plot(sd.index, sd.values, "o-", ms=3.5, color=col,
                     label=f"TL={tl:g}")
    axes[1].plot(sd.index, sd.values[0] * np.sqrt(sd.index[0] / sd.index),
                 color=C["grey"], ls="--", lw=1.2, label=r"$\propto L^{-1/2}$")
    axes[0].set_xscale("log"); axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[0].set_xlabel("sequence length $L$")
    axes[0].set_ylabel("median $I_{ss}$  (band: IQR)")
    title(axes[0], "A. Location vs $L$")
    axes[1].set_xlabel("sequence length $L$"); axes[1].set_ylabel("SD of $I_{ss}$")
    title(axes[1], r"B. Scale vs $L$ (dashed: $L^{-1/2}$)")
    axes[2].boxplot([sub[np.isclose(sub.TL, tl)].iss_A.values for tl in tls],
                    positions=tls, widths=[0.12 * t for t in tls],
                    showfliers=False, patch_artist=True, medianprops=dict(color="w"),
                    boxprops=dict(facecolor=C["blue"], alpha=0.55, lw=0.7),
                    whiskerprops=dict(color=C["grey"]),
                    capprops=dict(color=C["grey"]))
    axes[2].set_xscale("log")
    axes[2].set_xticks(tls); axes[2].set_xticklabels([f"{t:g}" for t in tls])
    axes[2].set_xlabel("total tree length"); axes[2].set_ylabel("$I_{ss}$")
    title(axes[2], "C. Location and spread vs tree length")
    for ax in axes[:2]:
        ax.legend(fontsize=6.5)
    save(fig, "fig1b_Lmap_profiles")


# ================================================================ Fig 2
def fig_iss_pdist():
    f = res(RES / "exp2_iss_pdist", "cells.csv")
    if not f.exists():
        return
    d = pd.read_csv(f)
    d["arm"] = d.cell.str.replace(r"_L\d+$", "", regex=True)
    d["fam"] = d.arm.str.split("_").str[0]
    sys.path.insert(0, str(ROOT / "src"))
    from issbench import exact, trees
    from issbench.models import build_Q

    Q = build_Q("JC69")
    PI = np.ones(4) / 4

    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2))
    colours = {4: C["black"], 8: C["blue"], 16: C["green"]}
    # invert p = 0.75 (1 - exp(-8 d/3)) for the per-branch length d, so the
    # exact curve is drawn against the same p-distance axis as the points
    pg = np.linspace(0.02, 0.7495, 120)
    br = -3.0 / 8.0 * np.log(np.clip(1 - 4 * pg / 3, 1e-12, None))
    for N, col in colours.items():
        s = d[d.fam == f"star{N}"]                    # shape-controlled: all stars
        if s.empty:
            continue
        for ax, conv in zip(axes[:2], ("A", "B")):
            ax.scatter(s.pbar, s.iss_A if conv == "A" else s.iss_B, s=5,
                       alpha=0.30, color=col, label=f"N={N}")
        cur = {conv: [] for conv in ("A", "B")}
        for d_br in br:
            node = trees.star(N, float(d_br))
            for conv in ("A", "B"):
                cur[conv].append(exact.E_Iss(N, node, Q, PI, convention=conv))
        for ax, conv in zip(axes[:2], ("A", "B")):
            ax.plot(pg, cur[conv], color=col, lw=1.6)
    # main arm: N=4 quartet with a 35% internal branch
    q4 = d[d.arm.str.startswith("q4")]
    axes[2].scatter(q4.pbar, q4.iss_A, s=5, alpha=0.30, color=C["black"],
                    label=r"$I_{ss}$ (A)")
    axes[2].scatter(q4.pbar, q4.iss_B, s=5, alpha=0.30, color=C["orange"],
                    label=r"$I_{ss}$ (B)")
    dgrid = np.exp(np.linspace(np.log(0.02), np.log(60.0), 90))
    pb, ia, ib = [], [], []
    for x in dgrid:
        r = 0.35 * x
        l = (x - r) / 4.0
        nd = trees.quartet(l, l, l, l, r)
        pb.append(0.75 * (1 - np.mean(np.exp(-4 * trees.pairwise_distances(nd)
                                             [np.triu_indices(4, 1)] / 3))))
        ia.append(exact.E_Iss(4, nd, Q, PI, convention="A"))
        ib.append(exact.E_Iss(4, nd, Q, PI, convention="B"))
    axes[2].plot(pb, ia, color=C["black"], lw=1.6, label="exact (A)")
    axes[2].plot(pb, ib, color=C["orange"], lw=1.6, label="exact (B)")
    for ax in axes:
        ax.set_xlabel("mean pairwise $p$-distance")
        ax.set_ylabel(r"$I_{ss}$")
        ax.set_xlim(0, 0.78)
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=6.8, loc="upper left")
    title(axes[0], "A. Star trees, convention A")
    title(axes[1], "B. Star trees, convention B")
    title(axes[2], r"C. Quartet, $f_{int}=0.35$ (N=4)")
    save(fig, "fig2_iss_vs_pdist")


# ================================================================ Fig 3
def fig_zone():
    rf = res(RES / "exp3_zone", "recovery.csv")
    cf = res(RES / "exp3_zone", "cells.csv")
    if not rf.exists():
        return
    rec = pd.read_csv(rf)
    wide = rec.pivot_table(index="t", columns="method", values="Ptrue")
    cells = pd.read_csv(cf)
    g = cells.groupby("cell").iss_A.agg(["mean", "std"]).reset_index()
    g["k"] = g["cell"].str.extract(r"t(\d+)").astype(int)
    ts = np.linspace(0.0, 1.0 - 0.05, 17)
    g["t"] = ts[g.k]
    g["tk"] = g["t"].round(6)
    wide.index = pd.Index(np.round(wide.index.to_numpy(), 6), name="t")
    g = g.sort_values("t")

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.5))
    mcol = {"ML": C["black"], "MPB": C["orange"], "NJ": C["blue"]}
    mk = {"ML": "o", "MPB": "^", "NJ": "D"}

    # A -- recovery
    for m, col in mcol.items():
        if m in wide.columns:
            axes[0].plot(wide.index, wide[m], mk[m] + "-", ms=4, lw=1.5,
                         color=col, label=m)
    axes[0].axhline(0.95, color=C["grey"], ls="--", lw=1)
    axes[0].set_xlabel("transfer parameter $t$  (Felsenstein to Farris)")
    axes[0].set_ylabel(r"$\Pr[\hat{T}=T]$")
    axes[0].legend(fontsize=7.5, loc="center left")
    title(axes[0], "A. Recovery along the transfer")

    # B -- Iss along the same path.  The EXPECTATION is flat to O(r); the
    # empirical means at 200 replicates scatter by more than that, so the exact
    # curve is overlaid rather than letting the reader read noise as structure.
    sys.path.insert(0, str(ROOT / "src"))
    from issbench import exact, trees as _t
    from issbench.models import build_Q as _bq
    P_, Q_, R_ = 1.0, 0.05, 0.02
    tgrid = np.linspace(0, P_ - Q_, 200)
    anal = np.array([exact.E_Iss(4, _t.quartet(P_ - x, Q_, P_, Q_ + x, R_),
                                 _bq("JC69"), np.ones(4) / 4, convention="A")
                     for x in tgrid])
    axes[1].plot(tgrid, anal, color=C["grey"], lw=2.4, zorder=1,
                 label=r"exact $\mathbb{E}[I_{ss}]$")
    axes[1].errorbar(g["t"], g["mean"], yerr=g["std"], fmt="o", ms=4.5,
                     color=C["purple"], capsize=2, zorder=3,
                     label="simulated mean $\\pm$ SD (200 reps)")
    gap = float(abs(anal[0] - anal[-1]))
    axes[1].set_xlabel("transfer parameter $t$")
    axes[1].set_ylabel(r"$I_{ss}$ (convention A)")
    lo = min(anal.min(), float((g["mean"] - g["std"]).min()))
    hi = max(anal.max(), float((g["mean"] + g["std"]).max()))
    pad = 0.12 * (hi - lo)
    axes[1].set_ylim(lo - pad, hi + pad)
    axes[1].legend(fontsize=7, loc="upper left")
    axes[1].annotate("exact endpoint gap %.4f $=O(r)$, $r=%.2f$" % (gap, R_),
                     xy=(0.03, 0.03), xycoords="axes fraction", fontsize=7.5,
                     color=C["grey"])
    title(axes[1], "B. The same path in $I_{ss}$ coordinates")

    # C -- the contrast, with the axis the reader needs for scale
    for m, col in mcol.items():
        if m not in wide.columns:
            continue
        mm = g.set_index("tk")["mean"]
        xs = [float(mm.loc[round(t, 6)]) for t in wide.index]
        axes[2].plot(xs, wide[m].values, "o", ms=4.5, color=col, alpha=0.8,
                     label=m)
    axes[2].axvline(0.6620, color=C["grey"], ls="--", lw=1)
    axes[2].text(0.6626, 18, "exact ceiling (conv. B)\n0.6620", fontsize=6.5,
                 color=C["grey"])
    axes[2].set_xlim(0.60, 0.70)
    axes[2].set_xlabel(r"$I_{ss}$ (convention A)")
    axes[2].set_ylabel(r"$\Pr[\hat{T}=T]$")
    axes[2].legend(fontsize=7.5, loc="center right")
    title(axes[2], "C. $\Pr[\hat{T}=T]$ against $I_{ss}$: one index value, "
                   "the full range of outcomes")
    save(fig, "fig3_zone_transfer")


# ================================================================ Fig 4
def fig_perturb():
    f = res(RES / "exp4_perturb", "cells.csv")
    if not f.exists():
        return
    d = pd.read_csv(f)
    d = d[d.L == 1000]

    def arm(prefix, num=True):
        s = d[d["cell"].str.startswith(prefix)]
        key = (s["cell"].str.extract(rf"{prefix}(inf|[\d.]+)")[0] if num
               else s["cell"].str.extract(rf"{prefix}([\d.]+)")[0])
        s = s.assign(x=key)
        g = s.groupby("x").iss_A.agg(["mean", "std"])
        return g

    fig, axes = plt.subplots(1, 5, figsize=(17.5, 3.5))

    # A -- alpha.  The homogeneous arm is tagged "alphainf", which a numeric
    # regex silently drops, so parse the suffix as text.
    a = arm("alpha")
    aord = ["inf", "4.0", "2.0", "1.0", "0.5", "0.25"]
    a = a.reindex([x for x in aord if x in a.index])
    axes[0].errorbar(range(len(a)), a["mean"], yerr=a["std"], fmt="o-", ms=4,
                     color=C["red"], capsize=2)
    axes[0].set_xticks(range(len(a)))
    axes[0].set_xticklabels(["$\infty$" if x == "inf" else x.rstrip("0").rstrip(".")
                             for x in a.index])
    axes[0].set_xlabel(r"$\alpha$"); axes[0].set_ylabel(r"$I_{ss}$")
    title(axes[0], r"A. $\Gamma$ shape $\alpha$")
    axes[0].annotate("homogeneous", xy=(0.02, 0.94), xycoords="axes fraction",
                     fontsize=6.5, color=C["grey"])

    for ax, (pre, lab, logx, col) in zip(axes[1:4], [
            ("kappa", "B. ts/tv $\kappa$", True, C["blue"]),
            ("fi", "C. internal fraction $f_i$", False, C["green"]),
            ("TL", "D. total tree length", True, C["purple"])]):
        s = arm(pre).reset_index()
        s["x"] = s["x"].astype(float)
        s = s.sort_values("x")
        ax.errorbar(s.x, s["mean"], yerr=s["std"], fmt="o-", ms=4, color=col,
                    capsize=2)
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(lab.split(". ", 1)[1]); ax.set_ylabel("$I_{ss}$")
        title(ax, lab)

    ax = axes[4]
    pis = [("pi_uniform", "uniform", C["black"]), ("pi_mild", "mild", C["orange"]),
           ("pi_atrich", "AT-rich", C["red"]), ("pi_extreme", "extreme", C["purple"])]
    vals, labs, cols = [], [], []
    for pre, lab, col in pis:
        s = d[d["cell"].str.startswith(pre)]
        if not s.empty:
            vals.append(s.iss_A.mean()); labs.append(lab); cols.append(col)
    ax.bar(range(len(vals)), vals, color=cols, width=0.6)
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=20, fontsize=7)
    ax.set_xlabel("composition"); ax.set_ylabel("$I_{ss}$")
    title(ax, "E. base composition")
    save(fig, "fig4_perturbation")


# ================================================================ Fig 5
def fig_condindep():
    """Is there anything left in Iss once pbar is controlled for?

    The test is not the correlation between the two statistics -- it is whether
    knowing Iss changes the predicted outcome at a FIXED pbar.  Each panel
    therefore bins on pbar, splits each bin by Iss, and plots the outcome for
    the low- and high-Iss halves of the bin.  A gap means Iss carries
    information beyond pbar; overlapping curves mean it does not.
    """
    d = RES / "exp6_condindep"
    f = res(d, "cells.csv")
    if not f.exists():
        return
    cells = pd.read_csv(f, usecols=["cell", "rep", "pbar", "iss_A"])
    targ = pd.read_csv(res(d, "targets.csv"), usecols=["cell", "rep", "pi_site"])
    df = targ.merge(cells, on=["cell", "rep"], how="inner")
    mpb = res(d, "recovery_targets.csv")
    ml = res(d, "recovery_ml.csv")
    if mpb.exists():
        df = df.merge(pd.read_csv(mpb)[["cell", "rep", "MPB"]],
                      on=["cell", "rep"], how="left")
    if ml.exists():
        df = df.merge(pd.read_csv(ml)[["cell", "rep", "ML"]],
                      on=["cell", "rep"], how="left")

    def by_iss_halves(y):
        """Bin on pbar, split each bin at its median Iss, return the two curves."""
        s = df.dropna(subset=[y]).copy()
        s["bin"] = pd.qcut(s.pbar, 30, labels=False, duplicates="drop")
        s["hi"] = s.iss_A > s.groupby("bin").iss_A.transform("median")
        g = s.groupby(["bin", "hi"]).agg(pbar=("pbar", "mean"), y=(y, "mean"),
                                         n=(y, "size")).reset_index()
        lo = g[~g.hi].sort_values("pbar")
        hi = g[g.hi].sort_values("pbar")
        return lo, hi

    panels = [("pi_site", "parsimony-informative fraction", C["red"]),
              ("MPB", "MPBoot recovery of the true quartet", C["orange"])]
    if "ML" in df.columns:
        panels.append(("ML", "ML (IQ-TREE) recovery of the true quartet",
                       C["blue"]))
    fig, axes = plt.subplots(1, len(panels) + 1,
                             figsize=(4.4 * (len(panels) + 1), 4.0))
    for ax, (y, lab, col) in zip(axes, panels):
        lo, hi = by_iss_halves(y)
        ax.plot(lo.pbar, lo.y, "o-", ms=3.5, lw=1.5, color=C["grey"],
                label="lower half of $I_{ss}$")
        ax.plot(hi.pbar, hi.y, "s-", ms=3.5, lw=1.5, color=col,
                label="upper half of $I_{ss}$")
        gap = float(np.nanmean(np.abs(lo.y.to_numpy() - hi.y.to_numpy())))
        ax.set_xlabel("mean pairwise $p$-distance")
        ax.set_ylabel(lab)
        ax.legend(fontsize=7, loc="best")
        title(ax, f"{'ABC'[panels.index((y, lab, col))]}. {lab}\n"
                  f"(mean gap {gap:.3f})")
    # final panel: what Iss ADDS, over the full range vs in the transition
    # band.  The full-range gains are ~0 for recovery, which is the number a
    # naive comparison reports; the band gains are where the decision lives.
    ax = axes[-1]
    cv = d / "condindep.csv"
    if cv.exists():
        c = pd.read_csv(cv)
        nm = {"MPB": "MPBoot", "ML": "ML (IQ-TREE)",
              "pi_site": "parsimony-info.", "pi_site_fraction": "parsimony-info."}
        c["lab"] = c.target.map(lambda t: nm.get(t, t))
        x = np.arange(len(c))
        ax.bar(x - 0.2, c.full_delta, 0.4, color=C["grey"],
               label=r"full range ($\bar p\in[0.1,0.75]$)")
        ax.bar(x + 0.2, c.band_delta, 0.4, color=C["green"],
               label=r"transition band ($\bar p\in[0.63,0.70]$)")
        ax.axhline(0, color=C["black"], lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(c.lab, fontsize=7.5)
        ax.set_ylabel(r"$\Delta$ AUC / $\Delta R^2$ from adding $I_{ss}$")
        ax.legend(fontsize=6.5, loc="upper left")
        title(ax, "D. What $I_{ss}$ adds, and where")
    save(fig, "fig5_condindep")


# ================================================================ Fig 6
def fig_empirical():
    """Where real data sit: the Xia & Lemey vertebrate COI partitions."""
    f = res(RES / "empirical", "summary.csv")
    if not f.exists():
        return
    d = pd.read_csv(f)
    d = d.set_index("partition").loc[["pos2", "pos1+2", "pos1", "all", "pos3"]]
    d = d.reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    x = np.arange(len(d))

    # A -- where each partition sits relative to the tabulated boundary
    ax = axes[0]
    ok = d.rf == 0
    ax.bar(x, d.iss_A, color=[C["green"] if o else C["red"] for o in ok],
           width=0.6)
    # The source reports BOTH boundaries for this data set at N = 8 and
    # L = 504: Iss.c = 0.7518 symmetrical, 0.6401 asymmetrical.
    for v, lab in ((0.7518, "symmetrical"), (0.6401, "asymmetrical")):
        ax.axhline(v, color=C["black"], ls="--", lw=1.2)
        ax.annotate(f"$I_{{ss.c}}$ = {v:.4f}, {lab}", xy=(1.05, v + 0.012),
                    fontsize=6.6, va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(d.partition, fontsize=8)
    ax.set_ylabel(r"$I_{ss}$ (convention A)")
    ax.set_ylim(0, 0.92)
    ax.set_xlabel("partition")
    title(ax, "A. The index, and the line it is read against")

    # B -- does the index track how much signal the partition carries?
    ax = axes[1]
    ax.scatter(d.n_informative, d.iss_A, s=70,
               color=[C["green"] if o else C["red"] for o in ok], zorder=3)
    for _, r in d.iterrows():
        ax.annotate(r.partition, (r.n_informative, r.iss_A),
                    textcoords="offset points", xytext=(7, 4), fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel("parsimony-informative sites")
    ax.set_ylabel(r"$I_{ss}$ (convention A)")
    title(ax, "B. Not a function of informativeness")

    # C -- recovery against the index
    ax = axes[2]
    # stagger the labels: the three recovered partitions sit on one line and
    # two of them are close in x
    dy = {"pos2": -12, "pos1+2": 14, "pos1": -12, "all": 10, "pos3": -12}
    for _, r in d.iterrows():
        ax.scatter(r.iss_A, r.rf, s=70,
                   color=C["green"] if r.rf == 0 else C["red"], zorder=3)
        ax.annotate(r.partition, (r.iss_A, r.rf), textcoords="offset points",
                    xytext=(8, dy[r.partition]), fontsize=7.5)
    ax.set_xlabel(r"$I_{ss}$ (convention A)")
    ax.set_ylabel("Robinson-Foulds distance to the reference tree")
    ax.set_ylim(-0.6, 3.2)
    title(ax, "C. Green: recovered; red: not")

    save(fig, "fig6_empirical")


# ================================================================ Fig 7
def fig_jarvis():
    """Where a genome-scale real alignment sits, whole and informative-only.

    Two subsets of the same concatenation: all complete columns, and the
    parsimony-informative columns within them, which is the restriction
    Duchene et al. recommend.  Same counts array, so both come from one pass.
    """
    f = res(RES / "empirical", "gigadb_windows.csv")
    fi = res(RES / "empirical", "gigadb_windows_pi.csv")
    g = res(RES / "empirical", "gigadb_global.csv")
    if not f.exists() or not g.exists():
        return
    d = pd.read_csv(f)
    di = pd.read_csv(fi) if fi.exists() else None
    glob = pd.read_csv(g)
    coi = res(RES / "empirical", "summary.csv")

    import sys as _s
    _s.path.insert(0, str(ROOT / "src"))
    from issbench import exact, trees
    from issbench.models import build_Q

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4))
    subsets = [("all complete columns", d, axes[0]),
               ("parsimony-informative sites", di, axes[1])]
    for title_txt, dd, row in subsets:
        if dd is None:
            continue
        tot = glob[glob.subset.str.contains("informative") ==
                   ("informative" in title_txt)].iloc[0]

        ax = row[0]
        ax.hist(dd.iss_A, bins=60, color=C["blue"], alpha=0.75,
                label=f"{len(dd):,} windows of 1000 sites")
        ax.axvline(float(tot.iss_A), color=C["red"], lw=1.8,
                   label=f"whole subset, $I_{{ss}}$={tot.iss_A:.3f}")
        ax.axvline(0.73, color=C["black"], ls="--", lw=1.4,
                   label="tabulated boundary ($N=8$)")
        if coi.exists():
            c = pd.read_csv(coi)
            for _, r in c.iterrows():
                ax.axvline(r.iss_A, color=C["orange"], lw=0.9, alpha=0.7)
            ax.annotate("vertebrate COI\ncodon partitions",
                        xy=(0.715, 0), xytext=(0.60, ax.get_ylim()[1] * 0.6),
                        fontsize=7, color=C["orange"])
        ax.set_xlabel("$I_{ss}$ per window")
        ax.set_ylabel("windows")
        ax.set_xlim(0, 0.85)
        ax.legend(fontsize=6.8, loc="upper right")
        panel(ax, "A" if "all" in title_txt else "C")
        ax.set_title(title_txt, fontsize=8.5, family="serif", pad=4)

        ax = row[1]
        ax.scatter(dd.pbar, dd.iss_A, s=3, alpha=0.25, color=C["grey"],
                   label=f"{len(dd):,} windows")
        ax.plot(dd.pbar.mean(), tot.iss_A, "*", ms=15, color=C["red"],
                label=f"whole subset ({tot.pbar:.3f}, {tot.iss_A:.3f})",
                zorder=5)
        pg = dd.pbar.quantile(np.linspace(0.01, 0.99, 40)).to_numpy()
        Q = build_Q("JC69"); N = int(tot.n_tax)
        ys = []
        for pp in pg:
            br = -3 / 8 * np.log(1 - 4 * pp / 3)
            ys.append(exact.E_Iss(N, trees.star(N, br), Q, np.ones(4) / 4,
                                  convention="A"))
        ax.plot(pg, ys, color=C["green"], lw=1.6,
                label=f"exact, JC69 star, $N$={N}")
        ax.set_xlabel("mean pairwise $p$-distance per window")
        ax.set_ylabel("$I_{ss}$ per window")
        ax.legend(fontsize=6.8, loc="upper left")
        panel(ax, "B" if "all" in title_txt else "D")

    save(fig, "fig7_jarvis")


def fig_lambda():
    """The pairwise lambda matrix localises what a global index cannot."""
    f = res(RES / "empirical", "lambda_matrix.csv")
    s = res(RES / "empirical", "lambda_matrix_summary.csv")
    if not f.exists():
        return
    m = pd.read_csv(f, index_col=0)
    d = pd.read_csv(s)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0),
                            gridspec_kw=dict(width_ratios=[1, 1.05]))
    ax = axes[0]
    im = ax.imshow(m.to_numpy(), cmap="YlOrRd", vmin=0, vmax=0.75)
    ax.set_xticks(range(len(m))); ax.set_xticklabels(m.columns)
    ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index)
    for a in range(len(m)):
        for b in range(len(m)):
            v = m.iloc[a, b]
            if not np.isnan(v):
                ax.text(b, a, f"{v:.2f}", ha="center", va="center",
                        fontsize=6.5, color="white" if v > 0.45 else "#333333")
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.046, label=r"$\lambda_{ij}$")
    panel(ax, "A")

    ax = axes[1]
    pairs = []
    for a in range(len(m)):
        for b in range(a + 1, len(m)):
            pairs.append(("with F" if "F" in (m.index[a], m.columns[b])
                          else "without F", m.iloc[a, b]))
    pdf = pd.DataFrame(pairs, columns=["group", "lam"])
    for k, (g, col) in enumerate((("without F", C["blue"]),
                                  ("with F", C["red"]))):
        v = pdf[pdf.group == g].lam
        ax.scatter(np.full(len(v), k) + np.random.uniform(-0.08, 0.08, len(v)),
                   v, s=28, color=col, alpha=0.85, label=g)
        ax.hlines(v.mean(), k - 0.25, k + 0.25, color=C["black"], lw=1.6)
        ax.annotate(f"mean {v.mean():.3f}", (k + 0.28, v.mean()), fontsize=7.5,
                    va="center")
    ax.axhline(float(d.lam_mean.mean()), color=C["grey"], ls="--", lw=1.2,
               label=f"mean over all pairs, {d.lam_mean.mean():.3f}")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["10 pairs\nwithout F", "5 pairs\nwith F"])
    ax.set_xlim(-0.5, 1.6)
    ax.set_ylabel(r"$\lambda_{ij}$")
    ax.legend(fontsize=7, loc="upper left")
    panel(ax, "B")
    save(fig, "fig8_lambda_matrix")


# ================================================================ Fig 9
def fig_tlc_alpha():
    """Critical tree length against rate heterogeneity, and the index at it.

    Section 2.6's claim in one picture.  Eleven values of alpha span two
    orders of magnitude, from no gamma at all to alpha = 0.5, with the
    no-gamma end sampled densely enough to show where the curve leaves it.
    The sweep stops at 0.5 because four gamma categories stop resolving a
    gamma below it (see the note in scripts/s14_tlc_alpha.py).
    """
    d = res(RES / "exp7_tlc_alpha", "pilot.csv")
    c = res(RES / "exp7_tlc_alpha", "critical.csv")
    if not d.exists() or not c.exists():
        return
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    p = pd.read_csv(d)
    cr = pd.read_csv(c).sort_values("alpha")
    # the same sweep on a Farris arrangement: identical tree length, identical
    # internal fraction, only the pendant split changed (scripts/s14b)
    c2 = RES / "exp7b_tlc_alpha_farris" / "critical.csv"
    cr2 = pd.read_csv(c2).sort_values("alpha") if c2.exists() else None

    # Eleven categories is past what a qualitative palette can carry, so the
    # ordering is a two-hue gradient and the scale a colour bar.  A colour sits
    # at its alpha's position on the ladder that panels B and C draw, with the
    # no-gamma case one doubling past the largest finite shape, so the same
    # alpha is in the same place in all three panels.  The earlier position,
    # log10(1 + 1/alpha), bunched everything above alpha = 8 into the last 3%
    # of the bar, which is exactly the end this sweep exists to resolve.
    XINF = 128.0

    def xa(a: float) -> float:
        return XINF if not np.isfinite(a) else float(a)

    alphas = sorted({float(x) for x in p.alpha.unique()})
    lo, hi = np.log10(xa(min(alphas))), np.log10(XINF)

    def pos(a: float) -> float:
        return (np.log10(xa(a)) - lo) / (hi - lo)

    cmap = LinearSegmentedColormap.from_list("alpha", [C["blue"], C["red"]])
    norm = Normalize(0.0, 1.0)
    col = {a: cmap(norm(pos(a))) for a in alphas}
    alab = {a: (r"$\alpha=\infty$" if not np.isfinite(a) else rf"$\alpha={a:g}$")
            for a in alphas}
    # only a handful of traces get a legend entry; the rest are read off the bar
    named = {np.inf, 4.0, 1.0, 0.5}

    def bar(ax):
        # every alpha on the bar would be eleven labels in an inch; the
        # intermediate values are read off the position, and the named traces
        # in the legend anchor the scale
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02)
        tick = [a for a in alphas if not np.isfinite(a) or a in (64.0, 16.0, 4.0, 1.0, 0.5)]
        cb.set_ticks([pos(a) for a in tick])
        cb.set_ticklabels([("∞" if not np.isfinite(a) else f"{a:g}")
                           for a in tick])
        cb.ax.tick_params(labelsize=6.5)
        cb.set_label(r"$\alpha$", fontsize=8)

    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2))
    ax = axes[0]
    # drawn from the light end up, so the wide curves at small alpha stay on top
    for a in alphas:
        sel = p[np.isclose(p.alpha, a)] if np.isfinite(a) else p[np.isinf(p.alpha)]
        g = sel.drop_duplicates(subset="tl").sort_values("tl")
        if g.empty:
            continue
        ax.plot(g.tl, g.ptrue, "o-", ms=2.6, lw=1.3, color=col[a],
                label=alab[a] if a in named else None)
    ax.axhline(0.95, color=C["grey"], ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("total tree length")
    ax.set_ylabel(r"$\Pr[\hat{T}=T]$  (ML)")
    ax.legend(fontsize=6.8, loc="lower left", title=r"$\Gamma$ shape",
              title_fontsize=7)
    panel(ax, "A")
    bar(ax)

    # infinity has no position on a log axis and matplotlib drops such points
    # without a word -- which silently removed the no-gamma reference from both
    # panels until xa() put it one doubling past the largest finite shape.
    ok = cr.dropna(subset=["tl_c"])

    def ladder(ax):
        ax.set_xscale("log")
        tick = [0.5, 1.0, 4.0, 16.0, 64.0, XINF]
        ax.set_xticks(tick)
        ax.set_xticklabels(["0.5", "1", "4", "16", "64", r"$\infty$"])
        ax.set_xlabel(r"$\Gamma$ shape $\alpha$")

    ax = axes[1]
    ax.plot([xa(a) for a in ok.alpha], ok.tl_c, ls="--", lw=0.9,
            color=C["grey"], zorder=1)
    for _, r in ok.iterrows():
        ax.plot(xa(r.alpha), r.tl_c, "o", ms=5.5, color=col[float(r.alpha)],
                zorder=3)
    if cr2 is not None:
        ok2 = cr2.dropna(subset=["tl_c"])
        ax.plot([xa(a) for a in ok2.alpha], ok2.tl_c, ls=":", lw=0.9,
                color=C["grey"], zorder=1)
        for _, r in ok2.iterrows():
            ax.plot(xa(r.alpha), r.tl_c, "s", ms=4.4, mfc="white",
                    mec=col[float(r.alpha)], mew=1.1, zorder=4)
    # the flat run from alpha = 8 to no gamma is one step wide; labelling each
    # point would collide, so it gets a single note
    for a in (4.0, 2.0, 1.0, 0.5):
        for _, r in ok[np.isclose(ok.alpha, a)].head(1).iterrows():
            ax.annotate(f"{r.tl_c:.0f}", (xa(r.alpha), r.tl_c), fontsize=7,
                        textcoords="offset points", xytext=(6, -3))
    flat = ok[ok.alpha >= 8]
    if len(flat):
        ax.annotate(f"$\\alpha\\geq8$:  {flat.tl_c.min():.1f}--"
                    f"{flat.tl_c.max():.1f}", (11.0, 22.0), fontsize=7.5,
                    color=C["black"], ha="center", va="center")
    ax.set_yscale("log")
    ax.set_ylabel("critical tree length")
    ladder(ax)
    panel(ax, "B")

    ax = axes[2]
    ax.plot([xa(a) for a in ok.alpha], ok.iss_c, ls="--", lw=0.9,
            color=C["grey"], zorder=1)
    for _, r in ok.iterrows():
        ax.plot(xa(r.alpha), r.iss_c, "o", ms=5.5, color=col[float(r.alpha)],
                zorder=3)
    for a, dx, dy in ((np.inf, -26, 3), (0.5, -6, 6), (1.0, 5, -11)):
        for _, r in ok[np.isclose(ok.alpha, a)].head(1).iterrows():
            ax.annotate(f"{r.iss_c:.3f}", (xa(r.alpha), r.iss_c), fontsize=7,
                        textcoords="offset points", xytext=(dx, dy))
    if cr2 is not None:
        ok2 = cr2.dropna(subset=["iss_c"])
        ax.plot([xa(a) for a in ok2.alpha], ok2.iss_c, ls=":", lw=0.9,
                color=C["grey"], zorder=1)
        for _, r in ok2.iterrows():
            ax.plot(xa(r.alpha), r.iss_c, "s", ms=4.4, mfc="white",
                    mec=col[float(r.alpha)], mew=1.1, zorder=4)
        sp1 = ok.iss_c.max() - ok.iss_c.min()
        sp2 = ok2.iss_c.max() - ok2.iss_c.min()
        ax.annotate(f"filled: four equal pendants, span {sp1:.3f}\n"
                    f"open: long pendants sisters, span {sp2:.3f}",
                    (0.03, 0.05), fontsize=7.0, color=C["black"],
                    xycoords="axes fraction", va="bottom")
    else:
        ax.annotate("span 0.0115, from $\\alpha=0.5$ to no $\\Gamma$",
                    (0.62, 0.45), fontsize=7.5, color=C["black"])
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"$I_{ss}$ at the critical tree length")
    ladder(ax)
    panel(ax, "C")
    save(fig, "fig9_tlc_alpha")


# ================================================================ Fig 10
def fig_theorem13():
    """Psi against p-bar: exact when the criterion holds, bent when it fails."""
    d = RES / "theorem13"
    f = d / "psi_affinity.csv"
    if not f.exists():
        return
    aff = pd.read_csv(f)
    files = {"HKY, uniform": "trace_HKY_uni.csv",
             "F81, skewed": "trace_F81_skew.csv",
             "HKY, skewed": "trace_HKY_skew.csv"}
    cols = {"HKY, uniform": C["blue"], "F81, skewed": C["green"],
            "HKY, skewed": C["red"]}

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    ax = axes[0]
    for name, fn in files.items():
        p = d / fn
        if not p.exists():
            continue
        t = pd.read_csv(p).sort_values("pbar")
        b, a = np.polyfit(t.pbar, t.psi, 1)
        ax.plot(t.pbar, t.psi, lw=2.4, color=cols[name], label=name)
        ax.plot(t.pbar, a + b * t.pbar, lw=0.9, ls="--", color=C["black"],
                alpha=0.75)
    ax.set_xlabel("mean pairwise $p$-distance")
    ax.set_ylabel(r"$\Psi=\langle\mathrm{tr}\,P\rangle-1$")
    ax.legend(fontsize=7.5, loc="upper right")
    panel(ax, "A")

    ax = axes[1]
    # residual in index units at N = 16: a residual in Psi moves the index by
    # (N-1) / (2 N h(pi)), and h depends on the composition, which differs
    # between the uniform and skewed cases
    import sys as _s
    _s.path.insert(0, str(ROOT / "src"))
    from issbench.entropy import h as _h
    hpi = {"HKY, uniform": _h(np.ones(4) / 4),
           "F81, skewed": _h(np.array([.4, .3, .2, .1])),
           "HKY, skewed": _h(np.array([.4, .3, .2, .1]))}
    for name, fn in files.items():
        p = d / fn
        if not p.exists():
            continue
        t = pd.read_csv(p).sort_values("pbar")
        b, a = np.polyfit(t.pbar, t.psi, 1)
        r = (t.psi - (a + b * t.pbar)) * (15.0 / (2 * 16 * hpi[name]))
        ax.plot(t.pbar, r, lw=2.0, color=cols[name], label=name)
    ax.axhline(0, color=C["grey"], lw=0.8)
    ax.set_xlabel("mean pairwise $p$-distance")
    ax.set_ylabel(r"residual, in $I_{ss}$ units ($N=16$)")
    ax.legend(fontsize=7.5, loc="upper left")
    panel(ax, "B")
    save(fig, "fig10_theorem13")


ALL = dict(fig_design=fig_design, fig_lmap=fig_lmap, fig_iss_pdist=fig_iss_pdist,
           fig_zone=fig_zone, fig_perturb=fig_perturb,
           fig_condindep=fig_condindep, fig_empirical=fig_empirical,
           fig_jarvis=fig_jarvis, fig_lambda=fig_lambda,
           fig_tlc_alpha=fig_tlc_alpha, fig_theorem13=fig_theorem13)

if __name__ == "__main__":
    style()
    want = [a for a in sys.argv[1:] if a in ALL] or list(ALL)
    for name in want:
        print(f"[{name}]")
        try:
            ALL[name]()
        except Exception as e:                                        # noqa: BLE001
            print(f"  FAILED: {type(e).__name__}: {e}")
