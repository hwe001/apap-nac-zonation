#!/usr/bin/env python3
"""Staggered-overdose recognition penalty: a virtual-trial sensitivity.

The main virtual trial draws presentation time independently of ingestion
pattern, which is optimistic for staggered patients: clinically, they present
late because a toxic *cumulative* dose is only recognised after repeated
ingestion. This script re-runs the standard-care arm for staggered patients
with presentation delayed by a penalty (default 12 h), keeping everything else
fixed (same seed -> same patients), and compares the pattern subgroup against
the unpenalised trial.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from virtual_trial import sample_patients, run_patient, ARMS  # noqa: E402
from nac_pk_calibration import cys_reference  # noqa: E402

RES = Path(__file__).resolve().parent.parent / "results"


def main(penalty_h=12.0, n_patients=300, seed=7):
    cys_ref, _, _ = cys_reference()
    pats = sample_patients(n_patients, seed)          # same seed -> same patients
    staggered = pats["staggered"].astype(bool)

    base = np.load(RES / "virtual_trial.npz")
    theta = {"assumed": float(base["theta_assumed"]),
             "measured": float(base["theta_measured"])}

    out = {}
    for scheme in ("assumed", "measured"):
        sev = np.empty(n_patients)
        for i in range(n_patients):
            pat = {k: v[i] for k, v in pats.items()}
            if staggered[i]:
                pat = dict(pat)
                pat["t_pres_h"] = pat["t_pres_h"] + penalty_h  # late recognition
            start = pat["t_pres_h"]                    # standard care
            cpk, _ = run_patient(scheme, pat, start, cys_ref)
            sev[i] = cpk > theta[scheme]
        out[scheme] = {
            "single": sev[~staggered].mean() * 100,
            "staggered": sev[staggered].mean() * 100,
            "all": sev.mean() * 100,
        }
        print(f"{scheme:9s} (+{penalty_h:.0f} h recognition penalty): "
              f"single {out[scheme]['single']:.1f}%  "
              f"staggered {out[scheme]['staggered']:.1f}%  "
              f"all {out[scheme]['all']:.1f}%")

    np.savez(RES / "vt_recognition.npz", penalty_h=penalty_h,
             **{f"sev_{s}_{k}": v for s in out for k, v in out[s].items()})
    return out


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 12.0)
