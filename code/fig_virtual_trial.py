#!/usr/bin/env python3
"""Fig 5 for paper 2: virtual-trial severity by presentation-time bin."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
FIGDIR = BASE / "figures"

ASSUMED_C = "#eb6834"
MEASURED_C = "#2a78d6"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e1e0d9"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titleweight": "bold", "axes.titlesize": 9.5, "figure.facecolor": "white",
    "axes.facecolor": "white", "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def main():
    d = np.load(RES / "virtual_trial.npz")
    t_pres = d["t_pres_h"]
    bins = [(0, 4), (4, 8), (8, 16), (16, 48)]
    labels = ["0–4 h", "4–8 h", "8–16 h", "> 16 h"]
    arms = [("none", "no NAC", "#8a8a8a"), ("std", "NAC at presentation", INK),
            ("late12", "NAC 12 h after presentation", "#955196")]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True)
    x = np.arange(len(bins))
    w = 0.26
    for ax, scheme, col in zip(axes, ("assumed", "measured"),
                               (ASSUMED_C, MEASURED_C)):
        theta = float(d[f"theta_{scheme}"])
        for j, (arm, lab, c) in enumerate(arms):
            vals = []
            for a, b in bins:
                m = (t_pres >= a) & (t_pres < b)
                vals.append((d[f"cpk_{scheme}_{arm}"][m] > theta).mean() * 100)
            ax.bar(x + (j - 1) * w, vals, w, color=c, label=lab,
                   edgecolor="white", lw=0.5)
            for xi, v in zip(x + (j - 1) * w, vals):
                ax.text(xi, v + 1, f"{v:.0f}", ha="center", fontsize=6.5,
                        color=INK2)
        arr = 42.3 - (d[f"cpk_{scheme}_std"] > theta).mean() * 100
        nnt = 100.0 / arr if arr > 0 else float("inf")
        ax.set_title(f"{scheme} gradients (NNT ≈ {nnt:.0f} for NAC at "
                     f"presentation)", fontsize=8.5, color=INK2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(True, lw=0.6, axis="y")
        ax.set_ylim(0, 80)
    axes[0].set_ylabel("severe-injury incidence (%)", fontsize=8)
    axes[0].legend(fontsize=7.5)
    fig.savefig(FIGDIR / "fig5_virtual_trial.png")
    plt.close(fig)
    print("wrote", FIGDIR / "fig5_virtual_trial.png")


if __name__ == "__main__":
    main()
