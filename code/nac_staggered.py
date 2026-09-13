#!/usr/bin/env python3
"""NAC treatment timing and staggered overdosing under assumed vs measured
zonation gradients.

Two extensions to the 16-hepatocyte model:

1. N-acetylcysteine (NAC): modeled as a temporary boost to hepatic GSH
   synthesis (cysteine delivery), bG -> bG + M * B_G, starting at t_start and
   lasting `dur` (default 21 h, the standard IV protocol length). The boost is
   uniform across cells (systemic cysteine supply) and M is expressed in units
   of basal synthesis capacity. We scan the treatment start time to map the
   predicted efficacy window under each gradient scheme.

2. Staggered ingestion: the instantaneous bolus is replaced by a first-order
   gut compartment (rate ka) with ingestion events added at specified times.
   This allows repeated/staggered overdose scenarios (Ghosh et al. 2026
   simulate single boluses only).

Metrics: peak pericentral protein adducts C_pk; NAC protection = 1 - C_pk(NAC)
/ C_pk(no NAC); pericentral GSH minimum and exhaustion time.

Units: time in days internally; doses in grams; NAC times in hours on the API
(converted internally).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
from intracellular_apap_model import (  # noqa: E402
    B_G, D_G, D_S, B_S, K_N, K_PSH, dose_to_p0,
)
from zonated_apap_model import build_parameters, K_CLEAR  # noqa: E402

FIGDIR = Path(__file__).resolve().parent.parent / "figures"
RES = Path(__file__).resolve().parent.parent / "results"

N = 16
# palette (same as make_figures)
ASSUMED_C = "#eb6834"
MEASURED_C = "#2a78d6"
INK2 = "#52514e"
GRID = "#e1e0d9"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.titleweight": "bold", "axes.titlesize": 9.5, "figure.facecolor": "white",
    "axes.facecolor": "white", "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def make_rhs_nac(p, nac=None, ka=0.0):
    """rhs for the 81-state system: [16 x (P,S,N,G,C), gut].

    nac = (t_start_days, dur_days, M): bG boosted by M*B_G during the window.
    ka > 0: gut absorption into all cells (P(0)=0, gut(0)=P0). ka = 0: bolus.
    """
    kS, kG, k450, kGSH, bG = p["kS"], p["kG"], p["k450"], p["kGSH"], p["bG"]

    if nac is None:
        t0, dur, boost = -1.0, 0.0, 0.0
    else:
        t0, dur, M = nac
        boost = M * B_G

    def rhs(t, y):
        cells = y[:80].reshape(N, 5)
        gut = y[80]
        P, S, Nq, G, C = cells[:, 0], cells[:, 1], cells[:, 2], cells[:, 3], cells[:, 4]
        bg_eff = bG + (boost if t0 <= t <= t0 + dur else 0.0)
        dP = -kS * S * P - kG * P - k450 * P + K_N * Nq + ka * gut
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P - K_N * Nq - kGSH * Nq * G - K_PSH * Nq
        dG = -kGSH * Nq * G + bg_eff - D_G * G
        dC = K_PSH * Nq - K_CLEAR * C
        return np.concatenate([np.stack([dP, dS, dN, dG, dC], axis=1).ravel(),
                               [-ka * gut]])

    return rhs


def integrate(scheme, dose_g, nac=None, ka=0.0, t_end=3.0, jumps=None, n_points=800):
    """Integrate; `jumps` = list of (t_days, gut_add_per_cell) for staggered dosing."""
    p = build_parameters(scheme)
    p0 = dose_to_p0(dose_g)
    g0 = p["bG"] / D_G
    y0 = np.concatenate([np.stack([np.full(N, p0 if ka == 0 else 0.0),
                                   np.full(N, B_S / D_S),
                                   np.zeros(N), g0, np.zeros(N)], axis=1).ravel(),
                         [p0 if ka > 0 else 0.0]])
    rhs = make_rhs_nac(p, nac=nac, ka=ka)
    t_eval = np.linspace(0, t_end, n_points)
    if not jumps:  # single segment
        sol = solve_ivp(rhs, (0, t_end), y0, method="Radau", t_eval=t_eval,
                        rtol=1e-5, atol=1e-8, max_step=0.02)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.t, sol.y[:80].reshape(N, 5, -1).transpose(2, 0, 1)

    # piecewise integration across ingestion jumps; drop each later segment's
    # first point (it duplicates the previous segment's boundary)
    seg_t, seg_Y = [], []
    y = y0.copy()
    bounds = [0.0] + sorted(t for t, _ in jumps) + [t_end]
    for a, b in zip(bounds[:-1], bounds[1:]):
        npts = max(2, int(n_points * (b - a) / t_end))
        sol = solve_ivp(rhs, (a, b), y, method="Radau",
                        t_eval=np.linspace(a, b, npts),
                        rtol=1e-5, atol=1e-8, max_step=0.02)
        if not sol.success:
            raise RuntimeError(sol.message)
        seg_t.append(sol.t)
        seg_Y.append(sol.y[:80])                 # (80, npts)
        y = sol.y[:, -1].copy()
        for tj, add in jumps:
            if abs(b - tj) < 1e-9:
                y[80] += add                     # gut jump at the ingestion time
    t = np.concatenate([seg_t[0]] + [s[1:] for s in seg_t[1:]])
    Yc = np.concatenate([seg_Y[0]] + [s[:, 1:] for s in seg_Y[1:]], axis=1)
    Y = Yc.reshape(N, 5, -1).transpose(2, 0, 1)  # (n_pts, 16, 5)
    return t, Y


def metrics(Y):
    Cpk_central = Y[:, -1, 4].max()
    Gmin_central = Y[:, -1, 3].min()
    g0c = None
    return {"Cpk": Cpk_central, "Gmin": Gmin_central}


# ---------------- Experiment 1: NAC efficacy window -------------------------
def nac_window():
    doses = (4.0, 16.0)
    starts_h = [0, 1, 2, 3, 4, 6, 8, 10, 12]
    dur_h, M = 21, 4
    out = {}
    for scheme in ("assumed", "measured"):
        out[scheme] = {}
        for dose in doses:
            _, Y0 = integrate(scheme, dose)
            base = metrics(Y0)["Cpk"]
            prot = []
            for ts in starts_h:
                _, Y = integrate(scheme, dose, nac=(ts / 24.0, dur_h / 24.0, M))
                prot.append((1 - metrics(Y)["Cpk"] / base) * 100)
            out[scheme][dose] = (base, prot)
            print(f"{scheme:9s} {dose:4.0f} g  Cpk_noNAC={base:.3e}  "
                  f"protection@0h={prot[0]:.0f}%  @8h={prot[6]:.0f}%  @12h={prot[-1]:.0f}%")
    RES.mkdir(exist_ok=True)
    np.savez(RES / "nac_window.npz",
             starts=starts_h,
             prot_assumed_4=out["assumed"][4.0][1], prot_assumed_16=out["assumed"][16.0][1],
             prot_measured_4=out["measured"][4.0][1], prot_measured_16=out["measured"][16.0][1])
    return out


def plot_nac_window(out):
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), sharey=True)
    starts_h = [0, 1, 2, 3, 4, 6, 8, 10, 12]
    for ax, dose in zip(axes, (4.0, 16.0)):
        for scheme, color in (("assumed", ASSUMED_C), ("measured", MEASURED_C)):
            prot = out[scheme][dose][1]
            ax.plot(starts_h, prot, "-o", color=color, lw=2, ms=5, label=scheme)
        ax.set_title(f"{dose:.0f} g", fontsize=9.5)
        ax.set_xlabel("NAC start time after ingestion (h)", fontsize=8)
        ax.set_xticks(starts_h)
        ax.tick_params(labelsize=8)
        ax.grid(True, lw=0.6)
        ax.axvspan(8, 12, color=GRID, alpha=0.5, zorder=0)
        ax.text(10, 2, "late\nwindow", fontsize=7, color=INK2, ha="center")
    axes[0].set_ylabel("reduction in peak pericentral\nadducts vs no NAC (%)", fontsize=8)
    axes[0].legend(fontsize=8, loc="upper right")
    fig.savefig(FIGDIR / "nac_window.png")
    plt.close(fig)


# ---------------- Experiment 2: staggered ingestion -------------------------
def staggered():
    ka = 60.0  # d^-1  (~2.5 h^-1, human APAP absorption)
    res, curves = {}, {}
    for scheme in ("assumed", "measured"):
        t1, Y1 = integrate(scheme, 16.0, ka=ka)             # single ingestion
        jumps1 = [(k * 2 / 24.0, dose_to_p0(4.0)) for k in range(1, 4)]
        t2, Y2 = integrate(scheme, 4.0, ka=ka, jumps=jumps1, t_end=3.0)   # 4 x 4 g q2h
        jumps2 = [(k * 1 / 24.0, dose_to_p0(2.0)) for k in range(1, 8)]
        t3, Y3 = integrate(scheme, 2.0, ka=ka, jumps=jumps2, t_end=3.0)   # 8 x 2 g q1h
        res[scheme] = {"single": metrics(Y1), "4x4g": metrics(Y2), "8x2g": metrics(Y3)}
        curves[scheme] = {"single": (t1, Y1), "4x4g": (t2, Y2), "8x2g": (t3, Y3)}
        print(scheme, {k: (f"Cpk={v['Cpk']:.3e}", f"Gmin={v['Gmin']:.2e}")
                       for k, v in res[scheme].items()})
    return res, curves


def plot_staggered(curves):
    """2x2: rows = scheme, cols = pericentral G(t) and C(t)."""
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 5.2), sharex=True)
    styles = {"single": ("-", 2.0), "4x4g": ("--", 1.8), "8x2g": (":", 1.8)}
    for r, scheme in enumerate(("assumed", "measured")):
        color = ASSUMED_C if scheme == "assumed" else MEASURED_C
        for c, (sp, unit) in enumerate(zip((3, 4), ("nmol", "nmol"))):
            ax = axes[r][c]
            for pat, (ls, lw) in styles.items():
                t, Y = curves[scheme][pat]
                ax.plot(t * 24, Y[:, -1, sp], ls=ls, lw=lw, color=color, label=pat)
            ax.set_ylabel(f"pericentral {'GSH' if sp == 3 else 'adducts'} ({unit})",
                          fontsize=8)
            ax.tick_params(labelsize=8)
            ax.grid(True, lw=0.6)
        axes[r][0].set_title(f"{scheme} gradients", fontsize=9, loc="left", color=INK2)
    for ax in axes[1]:
        ax.set_xlabel("time after first ingestion (h)", fontsize=8)
    axes[0][0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "staggered_patterns.png")
    plt.close(fig)


def main():
    FIGDIR.mkdir(exist_ok=True)
    print("=== NAC efficacy window (dur 21 h, M=4 x basal synthesis) ===")
    out = nac_window()
    plot_nac_window(out)
    print("\n=== staggered 16 g ingestion (ka = 60/d) ===")
    res, curves = staggered()
    plot_staggered(curves)
    print(f"\nfigures -> {(FIGDIR / 'nac_window.png').resolve()}")
    print(f"          {(FIGDIR / 'staggered_patterns.png').resolve()}")


if __name__ == "__main__":
    main()
