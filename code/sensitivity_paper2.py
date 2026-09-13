#!/usr/bin/env python3
"""Sensitivity and robustness sweeps for paper 2.

1. NAC boost magnitude M in {1, 2, 4} x basal GSH synthesis (the 21-h window
   length is held fixed): does the efficacy window depend on the (uncalibrated)
   boost size?
2. Absorption rate ka in {17, 34, 60} d^-1 (absorption half-lives ~1 h, 30 min,
   ~17 min): does the staggered-vs-single comparison depend on absorption?
3. Saturation robustness: replace the mass-action glucuronidation and oxidation
   terms with a saturable form v = k*P / (1 + P/Km), Km = P0(16 g)/r for
   r in {1, 3} (r -> inf is the linear model). Does the NAC window, and the
   assumed-vs-measured gap, survive saturation? (Ghosh et al. 2026 name
   Michaelis-Menten kinetics as their own next step.)
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
P0_16 = dose_to_p0(16.0)

ASSUMED_C = "#eb6834"
MEASURED_C = "#2a78d6"
INK2 = "#52514e"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "xtick.color": "#898781", "ytick.color": "#898781",
    "axes.titleweight": "bold", "axes.titlesize": 9.5, "figure.facecolor": "white",
    "axes.facecolor": "white", "grid.color": "#e1e0d9", "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight",
})


def make_rhs(p, nac=None, ka=0.0, sat_r=None):
    """rhs for [16 x (P,S,N,G,C), gut]. sat_r: Km = P0_16/sat_r for kG, k450."""
    kS, kG, k450, kGSH, bG = p["kS"], p["kG"], p["k450"], p["kGSH"], p["bG"]
    if nac is None:
        t0, dur, boost = -1.0, 0.0, 0.0
    else:
        t0, dur, M = nac
        boost = M * B_G
    km = P0_16 / sat_r if sat_r else None

    def rhs(t, y):
        cells = y[:80].reshape(N, 5)
        gut = y[80]
        P, S, Nq, G, C = cells[:, 0], cells[:, 1], cells[:, 2], cells[:, 3], cells[:, 4]
        bg_eff = bG + (boost if t0 <= t <= t0 + dur else 0.0)
        sat = 1.0 / (1.0 + P / km) if km else 1.0
        dP = -kS * S * P - kG * P * sat - k450 * P * sat + K_N * Nq + ka * gut
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P * sat - K_N * Nq - kGSH * Nq * G - K_PSH * Nq
        dG = -kGSH * Nq * G + bg_eff - D_G * G
        dC = K_PSH * Nq - K_CLEAR * C
        return np.concatenate([np.stack([dP, dS, dN, dG, dC], axis=1).ravel(),
                               [-ka * gut]])

    return rhs


def integrate(scheme, dose_g, nac=None, ka=0.0, sat_r=None, t_end=3.0,
              jumps=None, n_points=600):
    p = build_parameters(scheme)
    p0 = dose_to_p0(dose_g)
    g0 = p["bG"] / D_G
    y0 = np.concatenate([np.stack([np.full(N, p0 if ka == 0 else 0.0),
                                   np.full(N, B_S / D_S),
                                   np.zeros(N), g0, np.zeros(N)], axis=1).ravel(),
                         [p0 if ka > 0 else 0.0]])
    rhs = make_rhs(p, nac=nac, ka=ka, sat_r=sat_r)
    if not jumps:
        sol = solve_ivp(rhs, (0, t_end), y0, method="Radau",
                        t_eval=np.linspace(0, t_end, n_points),
                        rtol=1e-5, atol=1e-8, max_step=0.02)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.t, sol.y[:80].reshape(N, 5, -1).transpose(2, 0, 1)

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
        seg_Y.append(sol.y[:80])
        y = sol.y[:, -1].copy()
        for tj, add in jumps:
            if abs(b - tj) < 1e-9:
                y[80] += add
    t = np.concatenate([seg_t[0]] + [s[1:] for s in seg_t[1:]])
    Yc = np.concatenate([seg_Y[0]] + [s[:, 1:] for s in seg_Y[1:]], axis=1)
    return t, Yc.reshape(N, 5, -1).transpose(2, 0, 1)


def cpk(Y):
    return Y[:, -1, 4].max()


# ---------------- 1. M sweep -------------------------------------------------
def sweep_M():
    doses = (4.0, 16.0)
    starts_h = [0, 2, 4, 6, 8, 10, 12]
    Ms = (1.0, 2.0, 4.0)
    out = {}
    for scheme in ("assumed", "measured"):
        out[scheme] = {}
        for dose in doses:
            _, Y0 = integrate(scheme, dose)
            base = cpk(Y0)
            out[scheme][dose] = {"base": base, "M": {}}
            for M in Ms:
                prot = []
                for ts in starts_h:
                    _, Y = integrate(scheme, dose, nac=(ts / 24.0, 21 / 24.0, M))
                    prot.append((1 - cpk(Y) / base) * 100)
                out[scheme][dose]["M"][M] = prot
            print(f"M-sweep {scheme:9s} {dose:3.0f} g: " +
                  "  ".join(f"M={M:g}: {out[scheme][dose]['M'][M][4]:.0f}%@4h" for M in Ms))
    np.savez(RES / "sweep_M.npz", starts=starts_h, Ms=Ms,
             **{f"prot_{s}_{int(d)}g_M{int(M)}": out[s][d]["M"][M]
                for s in out for d in doses for M in Ms})
    return out, starts_h, Ms


def plot_M(out, starts_h, Ms):
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 5.6), sharex=True, sharey=True)
    blues = {1.0: "#b7d3f6", 2.0: "#5598e7", 4.0: "#1c5cab"}  # sequential ramp
    for i, scheme in enumerate(("assumed", "measured")):
        color = ASSUMED_C if scheme == "assumed" else None
        for j, dose in enumerate((4.0, 16.0)):
            ax = axes[i][j]
            for M in Ms:
                lw = 2.2 if M == 4.0 else 1.6
                c = color if scheme == "assumed" else blues[M]
                ax.plot(starts_h, out[scheme][dose]["M"][M], "-o", color=c,
                        lw=lw, ms=4.5, label=f"M = {M:g}")
            ax.set_title(f"{scheme}, {dose:.0f} g", fontsize=9, loc="left", color=INK2)
            ax.set_xticks(starts_h)
            ax.tick_params(labelsize=8)
            ax.grid(True, lw=0.6)
    for ax in axes[1]:
        ax.set_xlabel("NAC start time after ingestion (h)", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("reduction in peak\npericentral adducts (%)", fontsize=8)
    axes[0][0].legend(fontsize=8)
    fig.savefig(FIGDIR / "sweep_M.png")
    plt.close(fig)


# ---------------- 2. ka sweep ------------------------------------------------
def sweep_ka():
    kas = (17.0, 34.0, 60.0)
    rows = []
    for scheme in ("assumed", "measured"):
        for ka in kas:
            _, Y1 = integrate(scheme, 16.0, ka=ka)
            jumps = [(k * 2 / 24.0, dose_to_p0(4.0)) for k in range(1, 4)]
            _, Y2 = integrate(scheme, 4.0, ka=ka, jumps=jumps)
            r = {"scheme": scheme, "ka": ka,
                 "C_single": cpk(Y1), "C_stag": cpk(Y2),
                 "Gmin_single": Y1[:, -1, 3].min(), "Gmin_stag": Y2[:, -1, 3].min()}
            rows.append(r)
            print(f"ka-sweep {scheme:9s} ka={ka:4.0f}/d: single C={r['C_single']:.3e} "
                  f"stag C={r['C_stag']:.3e} (delta {100*(r['C_stag']/r['C_single']-1):+.1f}%)")
    return rows


# ---------------- 3. saturation robustness -----------------------------------
def sweep_sat():
    doses = (4.0, 16.0)
    starts_h = [0, 4, 8, 12]
    sat_rs = (None, 1.0, 3.0)
    out = {}
    for scheme in ("assumed", "measured"):
        out[scheme] = {}
        for dose in doses:
            out[scheme][dose] = {}
            for r in sat_rs:
                _, Y0 = integrate(scheme, dose, sat_r=r)
                base = cpk(Y0)
                prot = []
                for ts in starts_h:
                    _, Y = integrate(scheme, dose, nac=(ts / 24.0, 21 / 24.0, 4.0),
                                     sat_r=r)
                    prot.append((1 - cpk(Y) / base) * 100)
                out[scheme][dose][r] = {"base": base, "prot": prot}
                print(f"sat-sweep {scheme:9s} {dose:3.0f} g r={str(r):5s}: base={base:.3e}  "
                      f"prot@0h={prot[0]:.0f}%  @4h={prot[1]:.0f}%  @8h={prot[2]:.0f}%")
    np.savez(RES / "sweep_sat.npz", starts=starts_h,
             **{f"prot_{s}_{int(d)}g_r{str(r)}": out[s][d][r]["prot"]
                for s in out for d in doses for r in sat_rs})
    return out


def main():
    FIGDIR.mkdir(exist_ok=True)
    print("=== 1. NAC boost magnitude sweep (M = 1, 2, 4) ===")
    outM, starts_h, Ms = sweep_M()
    plot_M(outM, starts_h, Ms)
    print("\n=== 2. absorption-rate sweep (ka = 17, 34, 60 /d) ===")
    sweep_ka()
    print("\n=== 3. saturation robustness (linear vs Km-ratio r = 1, 3) ===")
    sweep_sat()
    print(f"\nfigures -> {(FIGDIR / 'sweep_M.png').resolve()}")


if __name__ == "__main__":
    main()
