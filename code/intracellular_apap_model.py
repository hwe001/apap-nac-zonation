#!/usr/bin/env python3
"""Intracellular APAP metabolism model (Reddyhoff et al. 2015), ported from
the original MATLAB implementation used for Means & Ho (2018/2019), 'A
spatial-temporal model for zonal hepatotoxicity of acetaminophen'.

Source MATLAB files (E:\\Google Drive\\Harvey Ho Paper\\2017\\APAP\\hepatotoxicity_code\\):
  - reddyhoff_apap.m          RHS function (ported directly below)
  - reddyhoff_apap_control.m  driver script (dosing, ICs, solver settings)

State variables (nmol): P = APAP, S = PAPS/sulphate conjugate, N = NAPQI,
G = glutathione (GSH), C = toxic protein adducts.

    dP/dt = -kS*S*P - kG*P - k450*P + kN*N
    dS/dt = -kS*S*P + bS - dS*S
    dN/dt = k450*P - kN*N - kGSH*N*G - kPSH*N
    dG/dt = -kGSH*N*G + bG - dG*G
    dC/dt = kPSH*N

Parameter values are taken from the MATLAB *code* (reddyhoff_apap.m), not
from Table 1 of the published paper -- cross-checking the two revealed the
published table's parameter *names* appear shuffled relative to their
values (most likely a transcription slip when the table was typed up for
submission). The code is internally consistent and reproduces the paper's
own stated dosing conversion (0.132 nmol for a 4 g dose) exactly, so it is
treated as ground truth here. See ../related_computational_models.md and
project conversation history for the cross-check.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

OUTPUT = Path("results")

# --- Parameters, exactly as in reddyhoff_apap.m ---------------------------
NMOL = 1e-12  # mol -> nmol rescaling factor

K_S = 2.26e14 * NMOL       # cell/(mol*day), rate PAPS converts APAP -> APAP-S
K_G = 2.99                 # 1/day, rate of APAP decay (glucuronidation)
K_450 = 0.315              # 1/day, cytochrome P450: APAP -> NAPQI
K_N = 0.0315                # 1/day, reverse rate NAPQI -> APAP
B_S = 2.65e-14 / NMOL       # mol/(cell*day), PAPS production by liver
D_S = 2.0                   # 1/day, PAPS decay
K_GSH = 1.6e18 * NMOL       # cell/(mol*day), NAPQI metabolised via GSH
K_PSH = 110.0                # 1/day, conjugation of NAPQI with protein
B_G = 1.374e-14 / NMOL       # mol/(cell*day), GSH production rate
D_G = 2.0                    # 1/day, GSH decay rate

# --- Dosing conversion, exactly as in reddyhoff_apap_control.m ------------
HEPATOCYTE_NUM = 1.6055e11
DOSAGE_PENETRATION = 0.8       # fraction of ingested dose reaching the liver
MOL_WEIGHT_APAP = 151.0        # g/mol


def dose_to_p0(dose_grams: float) -> float:
    """Ingested APAP dose (g) -> initial per-hepatocyte APAP concentration (nmol)."""
    p0 = (dose_grams * DOSAGE_PENETRATION / MOL_WEIGHT_APAP) / HEPATOCYTE_NUM
    return p0 / NMOL


@dataclass
class Solution:
    dose_grams: float
    t: np.ndarray
    P: np.ndarray
    S: np.ndarray
    N: np.ndarray
    G: np.ndarray
    C: np.ndarray


def rhs(t: float, y: np.ndarray) -> np.ndarray:
    P, S, N, G, C = y
    dP = -K_S * S * P - K_G * P - K_450 * P + K_N * N
    dS = -K_S * S * P + B_S - D_S * S
    dN = K_450 * P - K_N * N - K_GSH * N * G - K_PSH * N
    dG = -K_GSH * N * G + B_G - D_G * G
    dC = K_PSH * N
    return np.array([dP, dS, dN, dG, dC])


def run(dose_grams: float, t_end_days: float = 5.0, n_points: int = 2000) -> Solution:
    p0 = dose_to_p0(dose_grams)
    s0 = B_S / D_S
    g0 = B_G / D_G
    y0 = np.array([p0, s0, 0.0, g0, 0.0])

    t_eval = np.linspace(0, t_end_days, n_points)
    sol = solve_ivp(
        rhs, (0, t_end_days), y0, method="Radau",  # stiff solver, matches MATLAB ode15s
        t_eval=t_eval, rtol=1e-10, atol=1e-14,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed for dose={dose_grams}g: {sol.message}")

    return Solution(dose_grams, sol.t, sol.y[0], sol.y[1], sol.y[2], sol.y[3], sol.y[4])


def plot_comparison(sol_4g: Solution, sol_16g: Solution, output: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(9.5, 9), sharex=True)
    labels = [("P", "APAP"), ("S", "PAPS/APAP-S"), ("N", "NAPQI"), ("G", "GSH"), ("C", "Toxic adducts")]
    fields = ["P", "S", "N", "G", "C"]
    for ax, field, (short, long_) in zip(axes.flat, fields, labels):
        ax.plot(sol_4g.t, getattr(sol_4g, field), label="4 g (therapeutic)", color="#2a78d6")
        ax.plot(sol_16g.t, getattr(sol_16g, field), label="16 g (overdose)", color="#d03b3b")
        ax.set_title(f"{short}: {long_}", fontsize=10, fontweight="bold")
        ax.set_ylabel("nmol")
        ax.grid(alpha=0.3)
    axes.flat[-1].axis("off")
    axes[0, 0].legend(fontsize=9, frameon=False)
    for ax in axes[-1, :]:
        ax.set_xlabel("Time (days)")
    fig.suptitle("Intracellular APAP metabolism model (Reddyhoff 2015 parameters), Python port", fontsize=12)
    fig.tight_layout()
    fig.savefig(output / "intracellular_model_4g_vs_16g.png", dpi=200)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    p0_4g = dose_to_p0(4.0)
    p0_16g = dose_to_p0(16.0)
    print(f"p0 (4 g dose)  = {p0_4g:.6f} nmol  [paper states 0.132 nmol -- {'MATCH' if abs(p0_4g-0.132)<1e-3 else 'MISMATCH'}]")
    print(f"p0 (16 g dose) = {p0_16g:.6f} nmol")

    sol_4g = run(4.0)
    sol_16g = run(16.0)

    print()
    print("--- Sanity checks against the paper's qualitative Fig. 4 description ---")
    c_4g_final, c_16g_final = sol_4g.C[-1], sol_16g.C[-1]
    print(f"Final toxic adduct C: 4g={c_4g_final:.4g} nmol, 16g={c_16g_final:.4g} nmol, "
          f"ratio={c_16g_final/max(c_4g_final,1e-30):.1f}x  [paper: 'about two orders higher']")

    g_4g_min, g_16g_min = sol_4g.G.min(), sol_16g.G.min()
    print(f"Minimum GSH (G): 4g={g_4g_min:.4g} nmol, 16g={g_16g_min:.4g} nmol  "
          f"[paper: '16g depletion much faster']")

    s_4g_final, s_16g_final = sol_4g.S[-1], sol_16g.S[-1]
    s0 = B_S / D_S
    print(f"Final PAPS/APAP-S (S): 4g={s_4g_final:.4g} nmol, 16g={s_16g_final:.4g} nmol, "
          f"baseline s0={s0:.4g} nmol  [paper: 'sulphation saturated at 16g, [S] does not increase']")

    plot_comparison(sol_4g, sol_16g, OUTPUT)
    print(f"\nPlot written to {(OUTPUT / 'intracellular_model_4g_vs_16g.png').resolve()}")


if __name__ == "__main__":
    main()
