"""Figures for Part 2.

House style is Part 1's, copied rather than re-invented: Arial axis labels,
serif panel letters, muted palette, thin grey spines, light grid.  A Part 2
figure placed next to a Part 1 figure should look like it belongs.

    python Part_2/scripts/s41_figures.py            # all
    python Part_2/scripts/s41_figures.py fig1       # one

Every figure reads a CSV the analysis scripts wrote; none of them computes
anything.  If a CSV is missing the figure is skipped with a note rather than
silently drawn empty.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent

RES = PART2 / "results"
FIGS = PART2 / "figures"
FIGS.mkdir(exist_ok=True)

C = dict(black="#1a1a1a", red="#B22222", orange="#E8871E", blue="#3B7EA1",
         green="#2E8B57", purple="#7B4F9D", grey="#8C8C8C")
CREAM = "#FDF3E0"


def style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
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
    ax.set_title(s, loc="left", fontsize=11, fontweight="bold",
                 family="serif", color="#1a1a1a", pad=5)


def save(fig, name):
    p = FIGS / f"{name}.png"
    try:
        fig.tight_layout()
    except Exception:                                                # noqa: BLE001
        pass
    fig.savefig(p, bbox_inches=None)
    plt.close(fig)
    print(f"  wrote {p.relative_to(ROOT)}")


def need(*paths):
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        for p in missing:
            print(f"  skip: {Path(p).name} not found")
        return False
    return True


# ================================================================ Fig 1
def fig1():
    """Iss is constant; power is not."""
    f = RES / "iss_null_path" / "iss_null_path.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 5.0), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1, 1.6]))

    ax1.axhspan(d.iss.min() - 2e-15, d.iss.max() + 2e-15,
                color=CREAM, zorder=0)
    ax1.axhline(d.iss.iloc[0], color=C["red"], lw=1.6, zorder=3)
    ax1.set_ylim(0.55, 0.72)
    ax1.set_ylabel(r"$I_{ss}$")
    ax1.annotate(f"{d.iss.iloc[0]:.8f}  (spread $3\\times10^{{-15}}$)",
                 xy=(0.5, d.iss.iloc[0]), xytext=(0.5, 0.575),
                 ha="center", fontsize=8, color=C["red"],
                 arrowprops=dict(arrowstyle="-", color=C["red"], lw=0.8))
    panel(ax1, "A")

    ax2.plot(d.f, d.p_correct, "o-", color=C["black"], ms=4, lw=1.4,
             label=r"$\Pr[\mathrm{correct}]$")
    ax2.plot(d.f, d.mc, "s", color=C["blue"], ms=3.4, mfc="none",
             label="resampling check")
    ax2.axhline(1 / 3, color=C["grey"], ls="--", lw=1.0)
    ax2.text(0.30, 0.345, "chance", fontsize=8, color=C["grey"])
    ax2.set_xlabel("$f$  (fraction of columns shuffled within columns)")
    ax2.set_ylabel(r"$\Pr[\mathrm{correct}]$")
    ax2.set_ylim(0.25, 1.06)
    ax2.legend(loc="lower right", fontsize=8)
    panel(ax2, "B")

    save(fig, "fig1_iss_flat_power_not")


# ================================================================ Fig 2
def fig2():
    """Same Iss, two genuine trees."""
    f = RES / "matched_iss" / "matched_iss.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(figsize=(5.2, 3.9))
    ax.plot(d.iss, d.z_farris, "o-", color=C["red"], ms=5, lw=1.5,
            label="Farris   $a=b=1.0,\\ c=d=0.05$")
    ax.plot(d.iss, d.z_balanced, "s-", color=C["blue"], ms=5, lw=1.5,
            label="balanced  $a=b=c=d=0.1$")
    for _, r in d.iterrows():
        ax.annotate("", xy=(r.iss, r.z_balanced), xytext=(r.iss, r.z_farris),
                    arrowprops=dict(arrowstyle="-", color=C["grey"],
                                    lw=0.6, ls=":"))
    ax.axhline(1.645, color=C["grey"], ls="--", lw=1.0)
    ax.text(0.175, 1.9, "95% recovery", fontsize=8, color=C["grey"])
    ax.set_xlabel(r"$I_{ss}$  (matched exactly between the two)")
    ax.set_ylabel(r"$z=\Delta\sqrt{L}/\sqrt{v_2}$")
    ax.set_yscale("log")
    ax.legend(loc="upper left", fontsize=8)
    panel(ax, "")
    save(fig, "fig2_matched_iss")


# ================================================================ Fig 3
def fig3():
    """The exponent is a ratio of two measurable rates."""
    f = RES / "factors" / "factors.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    d = d[np.isfinite(d.p_fit) & (d.p_fit > 0)]
    fig = plt.figure(figsize=(7.6, 3.6))
    gs = fig.add_gridspec(1, 3, wspace=0.46, bottom=0.22, top=0.86,
                          left=0.07, right=0.98)

    # (a) predicted vs fitted
    ax = fig.add_subplot(gs[0, 0])
    lim = [0.03, 0.8]
    ax.plot(lim, lim, color=C["grey"], lw=0.9, ls="--")
    ax.plot(d.p_pred, d.p_fit, "o", color=C["black"], ms=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"$p$ predicted")
    ax.set_ylabel(r"$p$ fitted")
    err = (d.p_pred - d.p_fit).abs() / d.p_fit
    ax.text(0.05, 0.93, f"mean |rel. err.| = {err.mean():.3f}",
            transform=ax.transAxes, fontsize=8, va="top")
    panel(ax, "A")

    # (b) the factor plane
    ax = fig.add_subplot(gs[0, 1])
    for pv, ls in ((0.1, ":"), (0.2, ":"), (0.4, ":")):
        xs = np.linspace(0.4, 5.4, 50)
        ax.plot(xs, pv * xs, ls, color=C["grey"], lw=0.8)
        ax.text(5.0, pv * 5.0, f"$p$={pv}", fontsize=7, color=C["grey"],
                va="bottom", ha="right")
    ax.plot(d.lam_Delta, d.lam_Iss, "o", color=C["black"], ms=5)
    ref = d.iloc[0]
    for _, r in d.iloc[1:].iterrows():
        ax.annotate("", xy=(r.lam_Delta, r.lam_Iss),
                    xytext=(ref.lam_Delta, ref.lam_Iss),
                    arrowprops=dict(arrowstyle="->", color=C["orange"],
                                    lw=1.0, shrinkA=4, shrinkB=4))
    ax.set_xlabel(r"$\lambda_\Delta$")
    ax.set_ylabel(r"$\lambda_{I_{ss}}$")
    panel(ax, "B")

    # (c) curvature
    f3 = RES / "factors" / "local_ratio.csv"
    ax = fig.add_subplot(gs[0, 2])
    if f3.exists():
        lr = pd.read_csv(f3)
        for i, (cfg, sub) in enumerate(lr.groupby("config")):
            col = [C["red"], C["blue"], C["green"]][i % 3]
            ax.plot(sub.p_local, sub.ratio, "o", color=col, ms=4,
                    label=cfg.split()[0])
        lo = lr[["p_local", "ratio"]].min().min() * 0.9
        hi = lr[["p_local", "ratio"]].max().max() * 1.1
        ax.plot([lo, hi], [lo, hi], color=C["grey"], lw=0.9, ls="--")
        ax.legend(fontsize=7, loc="upper left")
    ax.set_xlabel(r"local $p$ from $L$")
    ax.set_ylabel(r"$\lambda_{I_{ss}}/\lambda_\Delta$ at $s^*(L)$")
    panel(ax, "C")

    save(fig, "fig3_exponent_decomposition")


# ================================================================ Fig 4
def fig4():
    """One quantity, three appearances.

    Panel C is the point of the figure: with the pendant sum, the internal
    branch and the total tree length all held fixed, moving only the
    arrangement drives Delta down by a factor of 28 in the Farris family and
    leaves it flat in the Felsenstein family -- which is what the Farris zone
    is.  Iss, shown alongside on the same axis for scale, moves by 2.6% over
    the same range.
    """
    f = RES / "farris_slowest" / "farris_slowest.csv"
    g = RES / "shared_origin" / "shared_origin.csv"
    if not need(f, g):
        return
    d = pd.read_csv(f)
    so = pd.read_csv(g)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.4, 3.2),
                                        sharex=True)
    cols = {"Farris": C["red"], "Felsenstein": C["blue"]}
    for fam, sub in d.groupby("family"):
        sub = sub.sort_values("r")
        ax1.plot(sub.r, sub.p, "o-", color=cols[fam], ms=5, lw=1.5, label=fam)
        ax2.plot(sub.r, sub.TL_crit, "s-", color=cols[fam], ms=5, lw=1.5,
                 label=fam)
    ax3.plot(so.r, so.delta_farris, "o-", color=C["red"], ms=5, lw=1.5,
             label="Farris")
    ax3.plot(so.r, so.delta_felsenstein, "s-", color=C["blue"], ms=5, lw=1.5,
             label="Felsenstein")
    ax3.set_yscale("log")
    ax3.annotate("$\\times28$", xy=(0.8, so.delta_farris.iloc[-1]),
                 xytext=(0.62, 3.5e-5), fontsize=8, color=C["red"],
                 arrowprops=dict(arrowstyle="->", color=C["red"], lw=0.9))
    ax3.set_xlabel("$r$  (pendant asymmetry; sum fixed)")
    ax3.set_ylabel(r"$\Delta=\mathrm{KL}(P\,\|\,\mathrm{star}^*)$")
    ax3.legend(fontsize=7, loc="lower left")

    ax1.set_xlabel("$r$")
    ax1.set_ylabel(r"$p$   (slower convergence $\downarrow$)")
    ax2.set_xlabel("$r$  (pendant asymmetry; sum fixed)")
    ax2.set_ylabel(r"$T_L^*$ at $L=3\times10^4$   (larger = easier)")
    panel(ax1, "A")
    panel(ax2, "B")
    panel(ax3, "C")
    save(fig, "fig4_arrangement")


# ================================================================ Fig 5
def fig5():
    """Where the criterion exists, and what slicing does to it."""
    f = RES / "blocking" / "blocking.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))

    # (a) threshold against block length, one curve per tree
    for i, (nm, sub) in enumerate(d.groupby("shape")):
        col = [C["red"], C["blue"]][i % 2]
        sub = sub.sort_values("L_block")
        ax1.plot(sub.L_block, sub.iss_c_block, "o-", color=col, ms=4, lw=1.4,
                 label=nm.split()[0] + " (block)")
        ax1.plot(sub.L_whole, sub.iss_c_whole, "^", color=col, ms=5,
                 mfc="none", label=nm.split()[0] + " (whole)")
        for _, r in sub.iterrows():
            ax1.annotate("", xy=(r.L_whole, r.iss_c_whole),
                         xytext=(r.L_block, r.iss_c_block),
                         arrowprops=dict(arrowstyle="->", color=col, lw=0.8,
                                         ls=":", alpha=0.8))
    ax1.set_xscale("log")
    ax1.set_xlabel("alignment length $L$")
    ax1.set_ylabel(r"$I_{ss.c}$")
    ax1.legend(fontsize=7, loc="lower right")
    panel(ax1, "A")

    # (b) how far the block mean sits above its own threshold
    r = d.sort_values("L_block")
    x = np.arange(len(r))
    ax2.bar(x - 0.18, r["ratio"], width=0.36, color=C["red"],
            label="shift / block sd")
    ax2.axhline(1, color=C["grey"], ls="--", lw=1.0)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{int(a)}\n{int(b)}" for a, b in
                         zip(r.L_block, r.L_whole)], fontsize=8)
    ax2.set_xlabel("$L_{block}$ / $L_{whole}$")
    ax2.set_ylabel("block mean above its\nown threshold (block sd)")
    ax2.legend(fontsize=7, loc="upper right")
    panel(ax2, "B")
    save(fig, "fig5_domain_and_slicing")


# ================================================================ Fig 6
def fig6():
    """Validation against resampling."""
    f = RES / "issc_definition" / "issc_definition.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    d["z"] = d["diff"] / d["iss_err"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.2))

    lab = {"a=b=1.0, c=d=0.05, e=0.02": "$e=0.02$",
           "a=b=1.0, c=d=0.05, e=0.10": "$e=0.10$",
           "a=b=1.0, c=d=0.20, e=0.05": "$c=d=0.20$",
           "a=b=0.6, c=d=0.05, e=0.02": "$a=b=0.6$"}
    for i, (nm, sub) in enumerate(d.groupby("shape")):
        col = [C["red"], C["blue"], C["green"], C["purple"]][i % 4]
        sub = sub.sort_values("L")
        ax1.errorbar(sub.L, sub["diff"], yerr=sub.iss_err, fmt="o",
                     color=col, ms=4, lw=1.0, capsize=2,
                     label=lab.get(nm, nm.split(",")[0]))
    ax1.axhline(0, color=C["black"], lw=0.9)
    ax1.set_xscale("log")
    ax1.set_xlabel("$L$")
    ax1.set_ylabel(r"$I_{ss.c}^{\rm MC}-I_{ss.c}^{\rm analytic}$")
    ax1.legend(fontsize=7, ncol=1)
    panel(ax1, "A")

    z = d.z.to_numpy()
    ax2.axvline(0, color=C["black"], lw=0.9)
    ax2.hist(z, bins=np.arange(-1.6, 1.7, 0.4), color=C["grey"],
             edgecolor="white", lw=0.6)
    ax2.axvline(z.mean(), color=C["red"], lw=1.6,
                label=f"mean {z.mean():+.3f}")
    ax2.axvline(np.median(z), color=C["blue"], lw=1.6, ls="--",
                label=f"median {np.median(z):+.3f}")
    ax2.set_xlabel("standardised residual")
    ax2.set_ylabel("count")
    ax2.legend(fontsize=7)
    npos = int((z > 0).sum())
    ax2.text(0.02, 0.96, f"sign test {npos}/{len(z)}  ($n<20$: not\n"
                         f"informative)", transform=ax2.transAxes,
             fontsize=7, va="top")
    panel(ax2, "B")
    save(fig, "fig6_validation")


ALL = {"fig1": fig1, "fig2": fig2, "fig3": fig3,
       "fig4": fig4, "fig5": fig5, "fig6": fig6}


def main() -> None:
    style()
    keys = sys.argv[1:] or list(ALL)
    for k in keys:
        if k not in ALL:
            print(f"  unknown figure {k}")
            continue
        ALL[k]()


if __name__ == "__main__":
    main()
