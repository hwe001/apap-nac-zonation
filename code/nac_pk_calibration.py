#!/usr/bin/env python3
"""Pharmacokinetics-anchored calibration of the NAC term.

Replaces the scenario parameter M (constant multiples of basal GSH synthesis)
with a protocol-driven input:

  PK   : plasma NAC from the standard 21-h IV protocol (150 mg/kg over 1 h,
         50 mg/kg over 4 h, 100 mg/kg over 16 h), first-order elimination with
         the reported IV half-life of ~2 h (Prescott et al.).
  PD   : a "cysteine-equivalent" pool fed by plasma NAC (pool t1/2 = 6 h,
         stated assumption, midpoint of the NAC t1/2 and hepatocyte GSH
         turnover), so that support persists a little beyond plasma levels.
  Boost: bG -> bG + M_max * B_G * (Cys / Cys_ref), i.e. the synthesis boost is
         proportional to the pool and peaks at M_max.

Calibration criterion (anchored to Lauterburg/Corcoran/Mitchell 1983, who found
NAC does NOT prevent initial GSH depletion but dramatically accelerates
recovery of detoxification capacity via de novo synthesis):

    at the 16 g dose with NAC started at 0 h, the pericentral hepatocyte's GSH
    recovers to >= 50% of its own pre-dose steady state by 12 h post-ingestion.

The smallest M_max meeting the criterion is the calibrated peak boost. The NAC
efficacy window is then re-run with the calibrated protocol-driven input and
compared with the constant-M scenario used previously.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
from intracellular_apap_model import (  # noqa: E402
    B_G, D_G, D_S, B_S, K_N, K_PSH, dose_to_p0,
)
from zonated_apap_model import build_parameters, K_CLEAR  # noqa: E402

N = 16
# --- PK parameters (literature) ---------------------------------------------
INFUS = ((0.0, 1 / 24, 150.0 * 24),      # (start_d, end_d, mg/kg/day)
         (1 / 24, 5 / 24, 50.0 * 24 / 4),
         (5 / 24, 21 / 24, 100.0 * 24 / 16))
KEL = np.log(2) / (2 / 24)               # plasma NAC t1/2 = 2 h (Prescott), 1/d
K_OUT = np.log(2) / (6 / 24)             # cysteine-equivalent pool t1/2 = 6 h
DOSE_TIMES_H = [0, 2, 4, 6, 8, 10, 12]
RES = Path(__file__).resolve().parent.parent / "results"


def infusion_rate(t):
    for a, b, r in INFUS:
        if a <= t < b:
            return r
    return 0.0


def cys_reference():
    """Peak of the cysteine-equivalent pool for the standard protocol."""
    def pk(t, y):
        A, Cys = y
        return [infusion_rate(t) - KEL * A, K_OUT * (A - Cys)]
    sol = solve_ivp(pk, (0, 21 / 24), [0, 0], method="LSODA",
                    t_eval=np.linspace(0, 21 / 24, 2000), rtol=1e-8, atol=1e-10)
    return sol.y[1].max(), sol.t, sol.y


def make_rhs_pk(p, M_max, cys_ref):
    kS, kG, k450, kGSH, bG = p["kS"], p["kG"], p["k450"], p["kGSH"], p["bG"]

    def rhs(t, y):
        cells = y[:80].reshape(N, 5)
        P, S, Nq, G, C = (cells[:, 0], cells[:, 1], cells[:, 2],
                          cells[:, 3], cells[:, 4])
        boost = M_max * B_G * (y[82] / cys_ref)
        bg_eff = bG + boost
        dP = -kS * S * P - kG * P - k450 * P + K_N * Nq
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P - K_N * Nq - kGSH * Nq * G - K_PSH * Nq
        dG = -kGSH * Nq * G + bg_eff - D_G * G
        dC = K_PSH * Nq - K_CLEAR * C
        return np.concatenate([
            np.stack([dP, dS, dN, dG, dC], axis=1).ravel(),
            [0.0,                                     # gut (unused: bolus mode)
             infusion_rate(t) - KEL * y[81],
             K_OUT * (y[81] - y[82])]])

    return rhs


def run_pk(scheme, dose_g, M_max, cys_ref, t_end=2.0, n_points=800):
    """16-cell model with bolus APAP and protocol-driven NAC (82 states)."""
    p = build_parameters(scheme)
    p0 = dose_to_p0(dose_g)
    g0 = p["bG"] / D_G
    y0 = np.concatenate([np.stack([np.full(N, p0), np.full(N, B_S / D_S),
                                   np.zeros(N), g0, np.zeros(N)], axis=1).ravel(),
                         [0.0, 0.0, 0.0]])
    sol = solve_ivp(make_rhs_pk(p, M_max, cys_ref), (0, t_end), y0,
                    method="Radau", t_eval=np.linspace(0, t_end, n_points),
                    rtol=1e-5, atol=1e-8, max_step=0.02)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y[:80].reshape(N, 5, -1).transpose(2, 0, 1)


def calibrate():
    cys_ref, t_pk, y_pk = cys_reference()
    print(f"PK: plasma-NAC peak {y_pk[0].max():.0f} mg/kg; "
          f"cysteine-pool peak {cys_ref:.0f} (pool t1/2 = 6 h)")

    baseline = B_G / D_G           # pericentral steady-state GSH, measured scheme
    target = 0.5 * baseline
    print(f"calibration criterion: pericentral G(12 h) >= {target:.2e} nmol "
          f"(50% of pre-dose steady state)\n")

    cal = None
    traj_t, traj_G, traj_M = None, {}, []
    for M_max in (0.0, 1.0, 2.0, 4.0, 8.0):     # 0 = no-NAC reference
        t, Y = run_pk("measured", 16.0, M_max, cys_ref, t_end=1.0)
        traj_t = t * 24.0                        # hours
        traj_G[M_max] = Y[:, -1, 3]
        traj_M.append(M_max)
        if M_max == 0.0:
            continue                             # reference: no criterion check
        g12 = np.interp(0.5, t, Y[:, -1, 3])   # pericentral G at t = 12 h
        ok = g12 >= target
        print(f"  M_max = {M_max:5.1f}: pericentral G(12 h) = {g12:.2e}  "
              f"({g12 / baseline * 100:5.1f}% of baseline)  "
              f"{'PASS' if ok else 'fail'}")
        if ok and cal is None:
            cal = M_max
    np.savez(RES / "nac_calibration.npz", t_h=traj_t, M=traj_M,
             baseline=baseline, **{f"G_M{int(M)}": traj_G[M] for M in traj_M})
    return cal, cys_ref


def window_with_pk(cal, cys_ref):
    """Protection vs start time under the calibrated protocol-driven NAC."""
    # start-time shift: NAC started later = skip the first portion of the
    # infusion; implemented by offsetting the infusion schedule via a delay.
    out = {}
    for scheme in ("assumed", "measured"):
        out[scheme] = {}
        for dose in (4.0, 16.0):
            _, Y0 = run_pk(scheme, dose, 0.0, cys_ref)  # M_max=0 -> no boost
            base = Y0[:, -1, 4].max()
            prot = []
            for ts in DOSE_TIMES_H:
                _, Y = run_pk_delayed(scheme, dose, cal, cys_ref, delay_h=ts)
                prot.append((1 - Y[:, -1, 4].max() / base) * 100)
            out[scheme][dose] = (base, prot)
            print(f"PK-NAC window {scheme:9s} {dose:3.0f} g: " +
                  "  ".join(f"{ts}h:{p:.0f}%" for ts, p in zip(DOSE_TIMES_H, prot)))
    np.savez(RES / "nac_window_pk.npz", starts=DOSE_TIMES_H,
             **{f"prot_{s}_{int(d)}g": out[s][d][1] for s in out for d in (4.0, 16.0)})
    return out


def run_pk_delayed(scheme, dose_g, M_max, cys_ref, delay_h):
    """Protocol-driven NAC with the infusion starting at `delay_h` post-dose."""
    p = build_parameters(scheme)
    kS, kG, k450, kGSH, bG = p["kS"], p["kG"], p["k450"], p["kGSH"], p["bG"]
    p0 = dose_to_p0(dose_g)
    g0 = p["bG"] / D_G
    delay = delay_h / 24.0

    def rhs(t, y):
        t_local = t - delay
        cells = y[:80].reshape(N, 5)
        P, S, Nq, G, C = (cells[:, 0], cells[:, 1], cells[:, 2],
                          cells[:, 3], cells[:, 4])
        boost = (M_max * B_G * (y[82] / cys_ref)) if t_local >= 0 else 0.0
        bg_eff = bG + boost
        I = infusion_rate(t_local) if t_local >= 0 else 0.0
        dP = -kS * S * P - kG * P - k450 * P + K_N * Nq
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P - K_N * Nq - kGSH * Nq * G - K_PSH * Nq
        dG = -kGSH * Nq * G + bg_eff - D_G * G
        dC = K_PSH * Nq - K_CLEAR * C
        return np.concatenate([
            np.stack([dP, dS, dN, dG, dC], axis=1).ravel(),
            [0.0, I - KEL * y[81], K_OUT * (y[81] - y[82])]])

    y0 = np.concatenate([np.stack([np.full(N, p0), np.full(N, B_S / D_S),
                                   np.zeros(N), g0, np.zeros(N)], axis=1).ravel(),
                         [0.0, 0.0, 0.0]])
    sol = solve_ivp(rhs, (0, 3.0), y0, method="Radau",
                    t_eval=np.linspace(0, 3.0, 800),
                    rtol=1e-5, atol=1e-8, max_step=0.02)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y[:80].reshape(N, 5, -1).transpose(2, 0, 1)


def main():
    cal, cys_ref = calibrate()
    print(f"\ncalibrated peak boost: M_max = {cal} x basal synthesis\n")
    if cal:
        window_with_pk(cal, cys_ref)


if __name__ == "__main__":
    main()
