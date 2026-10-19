"""Figures for the genome-scale section.

House style is s41's, so a figure added later still sits next to the others.

    python Part_2/scripts/s52_figures_genome.py            # all
    python Part_2/scripts/s52_figures_genome.py fig7       # one

Both figures read a CSV and compute nothing, as with s41.
"""
from __future__ import annotations

import importlib.util
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

_spec = importlib.util.spec_from_file_location("s41", str(HERE / "s41_figures.py"))
s41 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s41)
style, panel, save, need, C = s41.style, s41.panel, s41.save, s41.need, s41.C

RES = PART2 / "results"

ARM = {"CLEAN": C["blue"], "MID": C["orange"], "HIGH": C["red"],
       "ALL": C["grey"]}
ARM_LABEL = {"CLEAN": "CLEAN  ($I_{ss}<0.44$)",
             "MID": "MID  (0.44 ≤ $I_{ss}$ < 0.67)",
             "HIGH": "HIGH  ($I_{ss}\\geq 0.67$)",
             "ALL": "integrated"}


def fig8():
    """The band experiment: recovery against concatenated length."""
    f = RES / "caster_bands" / "bands_curve.csv"
    if not need(f):
        return
    d = pd.read_csv(f)
    style()
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.7), sharey=True)
    for ax, meth, letter in zip(
            axes, ["ML", "CASTER_site", "CASTER_pair"], "ABC"):
        for arm in ["CLEAN", "MID", "HIGH", "ALL"]:
            s = d[d.condition == arm].sort_values("L")
            if not len(s):
                continue
            ax.plot(s.L, s[meth], "o-" if arm != "ALL" else "s--",
                    color=ARM[arm], ms=3.4, lw=1.3,
                    mfc="white" if arm == "ALL" else ARM[arm],
                    label=ARM_LABEL[arm], zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(400, 1.4e6)
        ax.set_ylim(-1, 26)
        ax.axhline(0, color="#bbbbbb", lw=0.8, zorder=1)
        ax.set_xticks([500, 5_000, 50_000, 500_000])
        ax.set_xticklabels(["500", "5k", "50k", "500k"])
        ax.set_xlabel("concatenated length (bp)")
        name = {"ML": "IQ-TREE, JC", "CASTER_site": "CASTER-site",
                "CASTER_pair": "CASTER-pair"}[meth]
        panel(ax, f"{letter}   {name}")
    axes[0].set_ylabel("RF distance to the true tree")
    axes[0].legend(loc="upper right", fontsize=6.4, handlelength=1.6,
                   borderpad=0.4, labelspacing=0.25)
    save(fig, "fig8_band_experiment")


def fig7():
    """The admissible three-base weight closes onto -1/2."""
    f = RES / "caster_weights" / "window.csv"
    if not need(f):
        return
    d = pd.read_csv(f).sort_values("r")
    d = d[(d.r >= 4.0) & (d.r <= 26.0) & np.isfinite(d.window_lo)]
    style()
    fig, ax = plt.subplots(figsize=(3.6, 3.0))

    ax.fill_between(d.r, d.window_lo, d.window_hi,
                    color=C["green"], alpha=0.16, lw=0, zorder=1)
    ax.plot(d.r, d.window_lo, color=C["red"], lw=1.6,
            label="Farris lower bound", zorder=3)
    ax.plot(d.r, d.window_hi, color=C["blue"], lw=1.6,
            label="Felsenstein upper bound", zorder=3)
    ax.axhline(-0.5, color=C["black"], ls="--", lw=1.0, zorder=2)
    ax.text(4.2, -0.47, r"$w_3=-1/2$", fontsize=8, color=C["black"])
    ax.axvline(12, color=C["grey"], ls=":", lw=1.1, zorder=2)
    ax.text(12.3, -0.86, "end of\nthe Fig. S1 scan", fontsize=6.6,
            color=C["grey"], va="bottom")

    ax.set_xscale("log")
    ax.set_xlim(4, 26)
    ax.set_ylim(-1.15, -0.02)
    ax.set_xticks([4, 5, 6, 8, 10, 12, 15, 20, 25])
    ax.set_xticklabels(["4", "5", "6", "8", "10", "12", "15", "20", "25"])
    ax.minorticks_off()
    ax.set_xlabel("branch-length ratio $r$")
    ax.set_ylabel("admissible $w_3$")
    ax.legend(loc="lower right", fontsize=6.8, handlelength=1.6,
              borderpad=0.4, labelspacing=0.25)
    save(fig, "fig7_w3_window")


def main() -> None:
    want = [a for a in sys.argv[1:]] or ["fig7", "fig8"]
    for w in want:
        if w == "fig7":
            fig7()
        elif w == "fig8":
            fig8()


if __name__ == "__main__":
    main()
