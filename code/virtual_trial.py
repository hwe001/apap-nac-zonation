#!/usr/bin/env python3
"""Virtual (in silico) clinical trial: NAC treatment of acetaminophen overdose
in a virtual population, under two zonation representations.

Design
------
Virtual patients are sampled from distributions of the key uncertain inputs:

  total dose          log-normal, median 14 g, 5th-95th ~6-30 g (clipped 2-40 g)
  presentation time   log-normal, median 5 h, heavy tail (clipped 0.5-48 h)
  ingestion pattern   40% staggered (2-8 doses spread over 2-24 h), else single
  enzyme capacity     log-normal multipliers on total (not spatial) activity:
                      CYP450 CV 35%, UGT 30%, SULT 25%, GST 25%, GSH synthesis
                      30%, absorption rate 30% (stated assumptions)

Three paired arms (the SAME virtual patients in every arm):

  no-NAC     placebo (possible only in silico)
  standard   NAC at presentation (protocol-driven, calibrated treatment)
  late12     NAC 12 h after presentation (missed/delayed treatment)

Two parallel "worlds": assumed (Gebhardt-derived) and measured (spatial
proteomics) zonation gradients.

Endpoint anchoring: "severe injury" = peak pericentral adducts above a
threshold set per representation so that the no-treatment arm reproduces the
historical untreated severe-hepatotoxicity rate of 42.3% (Smilkstein et al.,
NEJM 1988; AST >= 1000 IU/L definition). Treatment arms are then compared
against that anchored threshold (absolute risk reduction, NNT, subgroup
effects). A secondary mechanistic endpoint is pericentral GSH exhaustion
(minimum < 10% of the patient's own pre-dose steady state).
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
from intracellular_apap_model import (  # noqa: E402
    B_G, D_G, D_S, B_S, K_N, K_PSH, dose_to_p0,
)
from zonated_apap_model import build_parameters, K_CLEAR  # noqa: E402
from nac_pk_calibration import KEL, K_OUT, cys_reference, infusion_rate  # noqa: E402

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
N = 16
M_MAX = 4.0                     # calibrated peak boost (main text)

# ---------------- population -------------------------------------------------
def sample_patients(n, seed=7):
    rng = np.random.default_rng(seed)
    dose = np.clip(np.exp(rng.normal(np.log(14.0), 0.4, n)), 2.0, 40.0)
    t_pres = np.clip(np.exp(rng.normal(np.log(5.0), 0.9, n)), 0.5, 48.0)
    staggered = rng.random(n) < 0.4
    n_doses = np.where(staggered, rng.integers(2, 9, n), 1)
    W = rng.uniform(2.0, 24.0, n)
    ln = lambda cv: np.exp(rng.normal(0.0, cv, n))
    return {
        "dose": dose, "t_pres_h": t_pres, "staggered": staggered,
        "n_doses": n_doses, "W": W,
        "m_cyp": ln(0.35), "m_ugt": ln(0.30), "m_sult": ln(0.25),
        "m_gst": ln(0.25), "m_bg": ln(0.30), "m_ka": ln(0.30),
    }


# ---------------- simulation -------------------------------------------------
def make_rhs(p, m_bg, nac_delay_d, M_max, cys_ref, ka):
    kS = p["kS"]; kG = p["kG"]; k450 = p["k450"]
    kGSH = p["kGSH"]; bG = p["bG"]
    treated = nac_delay_d is not None

    def rhs_full(t, y):
        cells = y[:80].reshape(N, 5)
        P, S, Nq, G, C = (cells[:, 0], cells[:, 1], cells[:, 2],
                          cells[:, 3], cells[:, 4])
        gut = y[80]
        boost = (M_max * B_G * m_bg * (y[82] / cys_ref)
                 if treated and t >= nac_delay_d else 0.0)
        I = infusion_rate(t - nac_delay_d) if treated and t >= nac_delay_d else 0.0
        dP = -kS * S * P - kG * P - k450 * P + K_N * Nq + ka * gut
        dS = -kS * S * P + B_S - D_S * S
        dN = k450 * P - K_N * Nq - kGSH * Nq * G - K_PSH * Nq
        dG = -kGSH * Nq * G + bG + boost - D_G * G
        dC = K_PSH * Nq - K_CLEAR * C
        return np.concatenate([
            np.stack([dP, dS, dN, dG, dC], axis=1).ravel(),
            [-ka * gut, I - KEL * y[81], K_OUT * (y[81] - y[82])]])

    return rhs_full


def run_patient(scheme, pat, nac_start_h, cys_ref, t_end=2.0):
    """One virtual patient, one arm. nac_start_h = absolute clock time (h since
    first ingestion) when the NAC infusion starts; None = no treatment.
    Returns (Cpk_pericentral, Gmin_frac)."""
    p = build_parameters(scheme)
    m_bg = pat["m_bg"]
    p = {k: v * m_bg for k, v in p.items()}          # total-capacity multiplier
    ka = 60.0 * pat["m_ka"]
    n_doses = int(pat["n_doses"])
    p0 = dose_to_p0(pat["dose"])
    per_dose = p0 / n_doses
    if n_doses == 1:
        jumps = None
        gut0 = p0
    else:
        W = pat["W"] / 24.0
        jumps = [(i * W / (n_doses - 1), per_dose) for i in range(1, n_doses)]
        gut0 = per_dose
    g0 = p["bG"] / D_G
    y0 = np.concatenate([np.stack([np.zeros(N), np.full(N, B_S / D_S),
                                   np.zeros(N), g0, np.zeros(N)], axis=1).ravel(),
                         [gut0, 0.0, 0.0]])
    delay = None if nac_start_h is None else nac_start_h / 24.0
    rhs = make_rhs(p, m_bg, delay, M_MAX, cys_ref, ka)

    def seg_solve(a, b, y, npts):
        sol = solve_ivp(rhs, (a, b), y, method="Radau",
                        t_eval=np.linspace(a, b, npts),
                        rtol=1e-5, atol=1e-8, max_step=0.02)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y

    if jumps is None:
        Y = seg_solve(0.0, t_end, y0, 300)[:80]
    else:
        bounds = [0.0] + sorted(t for t, _ in jumps) + [t_end]
        segs, y = [], y0.copy()
        for a, b in zip(bounds[:-1], bounds[1:]):
            npts = max(2, int(300 * (b - a) / t_end))
            Yb = seg_solve(a, b, y, npts)
            segs.append(Yb[:80])
            y = Yb[:, -1].copy()
            for tj, add in jumps:
                if abs(b - tj) < 1e-9:
                    y[80] += add
        Yc = np.concatenate([segs[0]] + [s[:, 1:] for s in segs[1:]], axis=1)
        Y = Yc
    Yc = Y.reshape(N, 5, -1).transpose(2, 0, 1)
    cpk = Yc[:, -1, 4].max()
    gmin_frac = Yc[:, -1, 3].min() / (p["bG"][-1] / D_G)
    return cpk, gmin_frac


# ---------------- trial driver -----------------------------------------------
ARMS = {"none": None, "std": 0.0, "late12": 12.0}


def main(n_patients=300, seed=7):
    RES.mkdir(exist_ok=True)
    cys_ref, _, _ = cys_reference()
    pats = sample_patients(n_patients, seed)
    t_start = time.time()

    results = {}
    for scheme in ("assumed", "measured"):
        results[scheme] = {}
        for arm in ARMS:
            cpk = np.empty(n_patients)
            gmin = np.empty(n_patients)
            t0 = time.time()
            for i in range(n_patients):
                pat = {k: v[i] for k, v in pats.items()}
                if arm == "none":
                    start = None                       # genuine placebo
                elif arm == "std":
                    start = pat["t_pres_h"]            # NAC at presentation
                else:                                  # late12
                    start = pat["t_pres_h"] + 12.0     # NAC 12 h after presentation
                cpk[i], gmin[i] = run_patient(scheme, pat, start, cys_ref)
            results[scheme][arm] = {"cpk": cpk, "gmin": gmin}
            print(f"{scheme:9s} {arm:7s}: {n_patients} pts in {time.time()-t0:.0f}s "
                  f"(median Cpk {np.median(cpk):.3e})", flush=True)

    # ---- anchoring and analysis -------------------------------------------
    summary = {}
    for scheme in ("assumed", "measured"):
        cpk_none = results[scheme]["none"]["cpk"]
        theta = np.quantile(cpk_none, 1.0 - 0.423)   # anchor: 42.3% untreated
        row = {"theta": theta}
        for arm in ARMS:
            sev = results[scheme][arm]["cpk"] > theta
            row[arm] = {"severity": sev.mean() * 100,
                        "gsh_exhaust": (results[scheme][arm]["gmin"] < 0.10).mean() * 100}
        arr = (row["none"]["severity"] - row["std"]["severity"]) / 100.0
        row["std"]["ARR_pp"] = row["none"]["severity"] - row["std"]["severity"]
        row["std"]["NNT"] = 1.0 / arr if arr > 0 else float("inf")
        arr_l = (row["none"]["severity"] - row["late12"]["severity"]) / 100.0
        row["late12"]["ARR_pp"] = row["none"]["severity"] - row["late12"]["severity"]
        row["late12"]["NNT"] = 1.0 / arr_l if arr_l > 0 else float("inf")
        summary[scheme] = row

    # ---- subgroup: severity vs presentation-time bin (standard arm) --------
    bins = [(0, 4), (4, 8), (8, 16), (16, 48)]
    subgroups = {}
    for scheme in ("assumed", "measured"):
        theta = summary[scheme]["theta"]
        t_pres = pats["t_pres_h"]
        sub = {}
        for arm in ARMS:
            sev = results[scheme][arm]["cpk"] > theta
            sub[arm] = [ (sev[(t_pres >= a) & (t_pres < b)]).mean() * 100
                         for a, b in bins ]
        sub["n_per_bin"] = [int(((t_pres >= a) & (t_pres < b)).sum()) for a, b in bins]
        subgroups[scheme] = sub

    # ---- pattern subgroup (standard arm) -----------------------------------
    pattern = {}
    for scheme in ("assumed", "measured"):
        theta = summary[scheme]["theta"]
        sev = results[scheme]["std"]["cpk"] > theta
        st = pats["staggered"].astype(bool)
        pattern[scheme] = {"single": sev[~st].mean() * 100,
                           "staggered": sev[st].mean() * 100,
                           "n_single": int((~st).sum()), "n_stag": int(st.sum())}

    # ---- save ---------------------------------------------------------------
    npz = {"seed": seed, "n": n_patients, "theta_assumed": summary["assumed"]["theta"],
           "theta_measured": summary["measured"]["theta"],
           "dose": pats["dose"], "t_pres_h": pats["t_pres_h"],
           "staggered": pats["staggered"]}
    for scheme in ("assumed", "measured"):
        for arm in ARMS:
            npz[f"cpk_{scheme}_{arm}"] = results[scheme][arm]["cpk"]
            npz[f"gmin_{scheme}_{arm}"] = results[scheme][arm]["gmin"]
    np.savez(RES / "virtual_trial.npz", **npz)

    with open(RES / "virtual_trial_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scheme", "arm", "severity_pct", "gsh_exhaustion_pct",
                    "ARR_pp_vs_none", "NNT"])
        for scheme in ("assumed", "measured"):
            for arm in ARMS:
                r = summary[scheme][arm]
                w.writerow([scheme, arm, f"{r['severity']:.1f}",
                            f"{r['gsh_exhaust']:.1f}",
                            f"{r.get('ARR_pp', ''):.1f}" if "ARR_pp" in r else "",
                            f"{r.get('NNT', ''):.1f}" if "NNT" in r else ""])

    print("\n=== anchored severity (theta set so no-NAC = 42.3%) ===")
    for scheme in ("assumed", "measured"):
        s = summary[scheme]
        print(f"{scheme:9s}: none {s['none']['severity']:.1f}%  "
              f"std {s['std']['severity']:.1f}% (ARR {s['std']['ARR_pp']:.1f} pp, "
              f"NNT {s['std']['NNT']:.1f})  "
              f"late12 {s['late12']['severity']:.1f}% (ARR {s['late12']['ARR_pp']:.1f} pp, "
              f"NNT {s['late12']['NNT']:.1f})  "
              f"| GSH exhaustion: none {s['none']['gsh_exhaust']:.0f}%  "
              f"std {s['std']['gsh_exhaust']:.0f}%")
    print("\n=== severity by presentation-time bin (standard arm) ===")
    labels = ["0-4h", "4-8h", "8-16h", "16h+"]
    for scheme in ("assumed", "measured"):
        sub = subgroups[scheme]
        print(f"{scheme:9s}: " + "  ".join(
            f"{lab}: none {sub['none'][i]:.0f}% std {sub['std'][i]:.0f}% "
            f"late {sub['late12'][i]:.0f}%"
            for i, lab in enumerate(labels)))
    print("\n=== pattern subgroup (standard arm) ===")
    for scheme in ("assumed", "measured"):
        p = pattern[scheme]
        print(f"{scheme:9s}: single {p['single']:.1f}% (n={p['n_single']})  "
              f"staggered {p['staggered']:.1f}% (n={p['n_stag']})")
    print(f"\ntotal wall time: {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    main(n)
