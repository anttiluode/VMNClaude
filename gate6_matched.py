"""Gate 6, post-hoc analysis (written AFTER seeing the pre-registered verdict).

The pre-registered matching target (single's answer error at a = 0.3, 0.0205) sits on the
memory's own error floor (~0.021, Gate 4), which the pair schemes approach but never reach
(best 0.0218). Matching at the floor is ill-posed. This script reports the damage each scheme
needs across a range of answer-error levels above the floor, by log-log interpolation along
each scheme's amplitude curve. Both damage measures: retained per-unit phase disturbance (rad)
and least-squares position shift.
"""
import json

import numpy as np

R = json.load(open("results/gate6_receipt.json"))
S, AMPS, SCH = R["summary_mean_sd"], R["setup"]["amps"], R["setup"]["schemes"]
TARGETS = (0.0225, 0.025, 0.03, 0.04, 0.06, 0.1)


def curve(sch, metric):
    return [(S[f"uniform|{sch}|{a:g}"]["answer_error"][0], S[f"uniform|{sch}|{a:g}"][metric][0]) for a in AMPS]


def cheapest(pts, target):
    """Least damage at which the scheme reaches answer error <= target, interpolating log-log between
    consecutive amplitudes. None if never reached."""
    best = None
    for (e1, d1), (e2, d2) in zip(pts, pts[1:]):
        for (e, d) in ((e1, d1), (e2, d2)):
            if e <= target:
                best = d if best is None else min(best, d)
        if (e1 - target) * (e2 - target) < 0:
            t = (np.log(target) - np.log(e1)) / (np.log(e2) - np.log(e1))
            d = float(np.exp(np.log(d1) + t * (np.log(d2) - np.log(d1))))
            best = d if best is None else min(best, d)
    return best


out = {}
for metric in ("phase_disturbance_rms", "position_shift"):
    out[metric] = {}
    for tgt in TARGETS:
        row = {s: cheapest(curve(s, metric), tgt) for s in SCH}
        row["single_over_sign_flip"] = row["single"] / row["sign_flip"] if row["single"] and row["sign_flip"] else None
        row["single_over_conjugate"] = row["single"] / row["conjugate"] if row["single"] and row["conjugate"] else None
        row["repeat_over_sign_flip"] = row["repeat"] / row["sign_flip"] if row["repeat"] and row["sign_flip"] else None
        out[metric][str(tgt)] = row
        print(metric[:14], tgt, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
json.dump(out, open("results/gate6_matched_posthoc.json", "w"), indent=2)
