#!/usr/bin/env python3
"""16-hepatocyte zonated APAP metabolism model.

Extends the single-cell Reddyhoff (2015) kinetics in `intracellular_apap_model.py`
across a periportal -> pericentral sinusoid of 16 hepatocytes, imposing zonal
gradients on the enzyme parameters. This is the computational core of the
"recalibration" paper: it lets us run the *assumed* (literature-derived)
gradient used by Means & Ho (2018) / LIU-Baylor (2026) side by side with the
*measured* gradient recalibrated from 2026 human spatial-omics.

State variables per hepatocyte i (nmol): P = APAP, S = PAPS/sulphate,
N = NAPQI, G = glutathione, C = toxic protein adducts. The 16 hepatocytes are
independent (no cell-to-cell transport -- matching Means & Ho 2018; transport
was added separately in Franiatte/Ho/Clarke 2019 and is out of scope here).

Zonated parameters and their measured directions (see
../measured_zonation_coefficients.md for the extracted numbers):

    parameter   assumed (Gebhardt/2018/2026)        measured (2026 human)
    --------    -------------------------------     ----------------------
    k450        pericentral (~5x)                   pericentral (confirmed, 5.4x)
    kG          pericentral (~9x)                   pericentral (confirmed, ~2.5x)
    kS          periportal  (~9x)                   UNZONATED  (flattened)
    kGSH        periportal  (~9x)                   PERICENTRAL (FLIPPED, 3.4x)
    bG (GSH)    periportal  (~9x)                   UNZONATED  (flattened)

The headline hypothesis: reversing kGSH (and flattening kS/bG) reduces the
predicted pericentral toxic-adduct hotspot -- i.e. the zone the model says is
most damaged is, in reality, the zone with the most GST detox capacity.
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
    K_S, K_G, K_450, K_N, B_S, D_S, K_GSH, K_PSH, B_G, D_G,
    dose_to_p0,
)

OUTPUT = Path("results")
N_HEPATOCYTES = 16

# Position along the sinusoid: 0.0 = periportal (hepatocyte #1),
# 1.0 = pericentral (hepatocyte #16).
POS = np.linspace(0.0, 1.0, N_HEPATOCYTES)

# First-order elimination of the adduct compartment (1/day). From the measured
# human APAP-protein-adduct elimination half-life t1/2 = 1.72 d
# (James et al. 2009, DMD 37:1779): k = ln2 / 1.72 = 0.403 d^-1; the paper's
# reported rate is 0.42 +/- 0.09 d^-1, used here. Not zonated -- adduct
# clearance is a whole-body process, not an enzyme gradient.
K_CLEAR = 0.42


def ramp(mult_portal: float, mult_central: float) -> np.ndarray:
    """Linear ramp of per-hepatocyte multipliers, portal -> central.

    Returns an array of length N_HEPATOCYTES. A multiplier > 1 means the enzyme
    is up-regulated relative to the uniform baseline; < 1 means down-regulated.
    """
    return mult_portal + (mult_central - mult_portal) * POS


def am_fold(F: float, central: bool = True) -> tuple:
    """Arithmetic-mean-preserving (mult_portal, mult_central) for a fold F.

    For a linear ramp the mean multiplier over the sinusoid is exactly
    (mult_portal + mult_central) / 2. To hold TOTAL enzyme abundance constant
    (equal to the uniform baseline) that mean must be 1, so the endpoints sum
    to 2. central=True: pericentral gene (central/portal = F); central=False:
    periportal gene (portal/central = F).
    """
    lo, hi = (2.0 / (1.0 + F), 2.0 * F / (1.0 + F)) if central else (2.0 * F / (1.0 + F), 2.0 / (1.0 + F))
    return (lo, hi)


# --- Zonal gradient schemes ---------------------------------------------
# Each entry maps a parameter to (mult_portal, mult_central).
#
# 'assumed': the qualitative Gebhardt (1992) gradient as operationalized by
#   Means & Ho 2018 and LIU-Baylor 2026 (their H1-H7, ~9x linear, CYP ~5x).
# 'assumed_kGSH_flipped': identical to 'assumed' except kGSH is reversed, to
#   isolate the single most important correction.
# 'measured': recalibrated from 2026 human spatial-omics.
SCHEMES = {
    "assumed": {
        # literature-derived (Gebhardt / LIU-Baylor H1-H7), ~9x linear, CYP ~5x;
        # arithmetic-mean-preserving so total enzyme == uniform baseline
        "k450": am_fold(5.0, True),    # CYP pericentral ~5x (H4)
        "kG":   am_fold(9.0, True),    # UGT pericentral ~9x (H3)
        "kS":   am_fold(9.0, False),   # SULT periportal ~9x (H2)
        "kGSH": am_fold(9.0, False),   # GST periportal ~9x (H6)  <- flips
        "bG":   am_fold(9.0, False),   # GSH synthesis periportal ~9x (H5)
    },
    "assumed_kGSH_flipped": {
        # identical to 'assumed' except kGSH reversed to pericentral (GSTA2 fold)
        "k450": am_fold(5.0, True),
        "kG":   am_fold(9.0, True),
        "kS":   am_fold(9.0, False),
        "kGSH": am_fold(1.16, True),   # FLIPPED (GSTA2 protein fold)
        "bG":   am_fold(9.0, False),
    },
    "measured": {
        # recalibrated from DVP binned intensity profiles (protein-level fold =
        # central-end / portal-end bins; NOT the robust-scaled coefficient).
        # Protein abundance is a proxy for catalytic capacity (assumption, stated).
        # k450 uses the isoform-abundance-weighted fold: CYP2E1/1A2/3A4 protein
        # folds 1.91/6.87/3.52 weighted 0.60/0.25/0.15 by their approximate
        # relative contribution to human APAP oxidation (Laine et al. 2009):
        #   0.60*1.91 + 0.25*6.87 + 0.15*3.52 = 3.39
        "k450": am_fold(3.39, True),   # weighted CYP (2E1 0.60, 1A2 0.25, 3A4 0.15)
        "kG":   am_fold(1.44, True),   # UGT2B7 (protein); UGT1A6 1.8x, UGT1A1 2.1x
        "kS":   (1.0, 1.0),            # SULT1A1 unzonated (fold 0.99x)
        "kGSH": am_fold(1.16, True),   # GSTA2 pericentral (protein fold 1.16x)
        "bG":   (1.0, 1.0),            # GCLC/GCLM/GSS unzonated (folds ~1.0x)
    },
}


def build_parameters(scheme):
    """Return per-hepatocyte parameter vectors for a scheme name or dict."""
    s = SCHEMES[scheme] if isinstance(scheme, str) else scheme
    return {
        "kS": K_S * ramp(*s["kS"]),
        "kG": K_G * ramp(*s["kG"]),
        "k450": K_450 * ramp(*s["k450"]),
        "kGSH": K_GSH * ramp(*s["kGSH"]),
        "bG": B_G * ramp(*s["bG"]),
    }


def make_rhs(p):
    kS, kG, k450, kGSH, bG = p["kS"], p["kG"], p["k450"], p["kGSH"], p["bG"]

    def rhs(t, y):
        y = y.reshape(N_HEPATOCYTES, 5)
        P = y[:, 0]
        S = y[:, 1]
        N = y[:, 2]
        G = y[:, 3]
        C = y[:, 4]
        dP = -kS * S * P - kG * P - k450 * P + K_N * N
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P - K_N * N - kGSH * N * G - K_PSH * N
        dG = -kGSH * N * G + bG - D_G * G
        dC = K_PSH * N - K_CLEAR * C
        return np.stack([dP, dS, dN, dG, dC], axis=1).ravel()

    return rhs


def run(scheme: str, dose_grams: float, t_end_days: float = 5.0, n_points: int = 1200,
        rtol: float = 1e-10, atol: float = 1e-14):
    """Integrate the 16-hepatocyte model and return (t, Y, params)."""
    p = build_parameters(scheme)
    p0 = dose_to_p0(dose_grams)
    s0 = B_S / D_S
    g0 = p["bG"] / D_G  # steady-state GSH consistent with each cell's bG
    y0 = np.stack([np.full(N_HEPATOCYTES, p0),
                   np.full(N_HEPATOCYTES, s0),
                   np.zeros(N_HEPATOCYTES),
                   g0,
                   np.zeros(N_HEPATOCYTES)], axis=1).ravel()

    t_eval = np.linspace(0, t_end_days, n_points)
    sol = solve_ivp(make_rhs(p), (0, t_end_days), y0, method="Radau",
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"integration failed for {scheme} {dose_grams}g: {sol.message}")
    Y = sol.y.reshape(N_HEPATOCYTES, 5, -1).transpose(2, 0, 1)  # (n_t, 16, 5)
    return sol.t, Y, p


def summarize(t, Y):
    """Final (t_end) adduct C profile and pericentral/periportal metrics."""
    C_final = Y[-1, :, 4]           # final toxic adducts per hepatocyte
    N_max = Y[:, :, 2].max(axis=0)  # peak NAPQI per hepatocyte
    portal = C_final[0]
    central = C_final[-1]
    ratio = central / portal if portal > 0 else float("inf")
    return C_final, N_max, {"portal_C": portal, "central_C": central, "ratio": ratio}


def perturb_assumed(param: str) -> dict:
    """'assumed' scheme with one parameter overridden by its 'measured' value."""
    s = dict(SCHEMES["assumed"])
    s[param] = SCHEMES["measured"][param]
    return s


def decompose(dose_grams: float) -> None:
    """One-at-a-time attribution: % change in pericentral C from moving a single
    parameter from its assumed to its measured gradient."""
    _, Y_base, _ = run("assumed", dose_grams)
    base_central = summarize(None, Y_base)[2]["central_C"]
    print(f"  One-at-a-time (baseline assumed central_C = {base_central:.4g}):")
    for param in ("k450", "kG", "kS", "kGSH", "bG"):
        _, Y_pert, _ = run(perturb_assumed(param), dose_grams)
        c = summarize(None, Y_pert)[2]["central_C"]
        print(f"    {param:6s} -> measured: central_C={c:.4g}  "
              f"({(c / base_central - 1) * 100:+.1f}%)")


def plot_profiles(results, dose_grams):
    """Spatial profile of final toxic adducts: assumed vs measured."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    for scheme, (t, Y, _) in results.items():
        C_final, _, _ = summarize(t, Y)
        ax.plot(np.arange(1, N_HEPATOCYTES + 1), C_final, "-o", ms=4, label=scheme)
    ax.set_xlabel("Hepatocyte index (1 = periportal, 16 = pericentral)")
    ax.set_ylabel("Final toxic protein adducts C (nmol)")
    ax.set_title(f"Final adduct profile, {dose_grams} g")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for scheme, (t, Y, _) in results.items():
        ax.plot(t, Y[:, -1, 4], label=scheme)  # pericentral cell #16 C(t)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Pericentral adducts C (nmol)")
    ax.set_title("Pericentral (#16) adduct time course")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

    fig.suptitle("Zonated APAP model: assumed vs measured gradients", fontsize=12)
    fig.tight_layout()
    fname = OUTPUT / f"zonated_{dose_grams:.0f}g_assumed_vs_measured.png"
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for dose in (4.0, 16.0):
        print(f"\n===== {dose:.0f} g dose =====")
        results = {}
        for scheme in ("assumed", "assumed_kGSH_flipped", "measured"):
            t, Y, _ = run(scheme, dose)
            results[scheme] = (t, Y, _)
            _, _, m = summarize(t, Y)
            print(f"  {scheme:22s}: central_C={m['central_C']:.4g}  "
                  f"portal_C={m['portal_C']:.4g}  central/portal ratio={m['ratio']:.2f}x")

        m_assumed = summarize(results["assumed"][0], results["assumed"][1])[2]
        m_flip = summarize(results["assumed_kGSH_flipped"][0],
                           results["assumed_kGSH_flipped"][1])[2]
        m_meas = summarize(results["measured"][0], results["measured"][1])[2]
        print(f"  -> kGSH flip alone changes pericentral C by "
              f"{(m_flip['central_C'] / m_assumed['central_C'] - 1) * 100:+.1f}%")
        print(f"  -> full recalibration changes pericentral C by "
              f"{(m_meas['central_C'] / m_assumed['central_C'] - 1) * 100:+.1f}%")

        decompose(dose)

        fname = plot_profiles(results, dose)
        print(f"  plot -> {fname.resolve()}")


if __name__ == "__main__":
    main()
