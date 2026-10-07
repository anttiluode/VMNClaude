"""Gate 5b - why antipodal sweeps stop cancelling while navigating (diagnostic for Gate 5).

Gate 5 Part B: with sweeps along the CURRENT heading, antipodal pairs cut read damage only
1.3x, against 35x when stationary. Diagnosis tested here: heading turns between the two
sweeps of a pair, so +D and -D' are not exact opposites. The fine scale (b = 8 pi) turns a
small offset mismatch s*dh into a large per-unit phase error.

Variant: the second sweep of each pair reuses the first sweep's heading ("paired").
Also reported: mean heading change between consecutive sweeps.
"""
import json

import numpy as np

import gate2_path as g2
import gate5_sweeps as g5

PATS = ("none", "same", "left_right", "antipodal", "left_right_paired", "antipodal_paired")


def run(sd, amp, every=2):
    _, d, b = g2.bank(sd)
    x, v = g2.trajectories(np.random.default_rng(300 + sd), g5.E_B)
    out = {}
    turn = []
    for pat in PATS:
        base = pat.replace("_paired", "")
        zp = np.sqrt(g5.MU) * np.exp(1j * b[:, None] * (d @ x[0].T))
        zu = zp.copy()
        nrng = np.random.default_rng(70 + sd)
        prng = np.random.default_rng(80 + sd)
        k, hprev = 0, None
        for t in range(g5.STEPS):
            lin = g5.MU + 1j * (b[:, None] * (d @ v[t].T) / (g5.SUB * g5.DT))
            if t % every == 0:
                h = v[t] / np.maximum(np.linalg.norm(v[t], axis=1, keepdims=True), 1e-12)
                if pat == "none" and hprev is not None:
                    turn.append(np.abs(np.angle((h[:, 0] + 1j * h[:, 1]) * np.conj(hprev[:, 0] + 1j * hprev[:, 1]))).mean())
                if pat.endswith("_paired") and k % 2 == 1:
                    h = hprev
                hprev = h
                zp = g5.ping(zp, g5.offsets(base, k, h, prng), amp, d, b)
                k += 1
            for _ in range(g5.SUB):
                ep, eu = g5.shared(g5.noise(nrng, zp.shape), zp, zu)
                zp, zu = g5.rk4(zp, lin, 0.0) + ep, g5.rk4(zu, lin, 0.0) + eu
        out[pat] = float(np.linalg.norm(g5.shift_of(zp, zu, d, b), axis=1).mean())
    out["mean_heading_change_rad_between_sweeps"] = float(np.mean(turn))
    return out


def main():
    res = {}
    for amp in (0.1, 0.3):
        per = [run(sd, amp) for sd in g5.SEEDS]
        res[f"a{amp:g}"] = {k: [float(np.mean([p[k] for p in per])), float(np.std([p[k] for p in per]))] for k in per[0]}
        print(amp, {k: round(v[0], 5) for k, v in res[f"a{amp:g}"].items()}, flush=True)
    json.dump(res, open("results/gate5b_receipt.json", "w"), indent=2)


if __name__ == "__main__":
    main()
