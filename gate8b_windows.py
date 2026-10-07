"""Gate 8b - diagnostic, written AFTER Gate 8's K8b failed: is the wire lossy, or the window?

Gate 8's one-wire frequency-multiplexed (FDM) listener reached answer error 0.124 even with no
measurement noise, against 0.022 for the full per-unit listener. Separating 40 frequencies on one
wire needs windows of at least N samples; at the bank's sample step that is 12 time units. The
hypothesis: the wire carries each unit's window-averaged state almost losslessly, and the loss
comes from the long window, which smears the ~1-time-unit amplitude transient carrying the
cos(alpha) half of the answer. Then a wire with more bandwidth (more samples per time unit)
could use short windows and recover the full listener's accuracy.

Test, with NO wire (each unit's state averaged directly over a window, i.e. a perfect
demodulator):
  avg_long    one 48-substep window before the ping, one after (Gate 8's timing)
  avg_short   one 4-substep window just before the ping; after it, windows over substeps 1-2,
              3-4, 5-8, 9-16. Equivalent to an FDM wire with 12x the bandwidth (40 samples per
              short window).
If avg_long ~ wire_fdm_sigma0 and avg_short ~ full, the diagnosis holds.
a = 0.2; 3 seeds; same banks, goals and noise streams as Gate 8.
"""
import json
import os

import numpy as np

import gate2_path as g2
import gate4_ping as g4
import gate6_restore as g6
import gate7_answer_undo as g7

MU = g2.MU
E, OFF, SEEDS, SPLIT_ID = g6.E, g6.OFF, g6.SEEDS, g7.SPLIT_ID
AMP = 0.2
L = 48
SHORT_POST = ((1, 2), (3, 4), (5, 8), (9, 16))


def feats_from(pre, posts):
    out = []
    for p in posts:
        out += [np.angle(p * np.conj(pre)), np.abs(p) - np.abs(pre)]
    return np.concatenate(out, 0).T


def run(sd, Z, q, rng):
    _, d, b = g2.bank(sd)
    phi = b[:, None] * (d @ q.T)
    zp, zu = Z.copy(), Z.copy()
    hist_pre = []
    for _ in range(L):
        hist_pre.append(zp)
        zp, zu = g6.step(zp, zu, rng)
    zp = zp + AMP * np.sqrt(MU) * np.exp(1j * phi)
    hist_post = []
    for _ in range(L):
        zp, zu = g6.step(zp, zu, rng)
        hist_post.append(zp)
    long_f = feats_from(np.mean(hist_pre, 0), [np.mean(hist_post, 0)])
    short_f = feats_from(np.mean(hist_pre[-4:], 0), [np.mean(hist_post[a - 1:b], 0) for a, b in SHORT_POST])
    return long_f, short_f


def main():
    per = {}
    for sd in SEEDS:
        data, Fl, Fs = {}, {}, {}
        for s in E:
            xx, v = g2.trajectories(np.random.default_rng(OFF[s] + sd), E[s])
            rng = np.random.default_rng(900 + OFF[s] + sd)
            ang = rng.uniform(0, 2 * np.pi, E[s]); rad = 0.35 * np.sqrt(rng.uniform(0, 1, E[s]))
            D = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
            data[s] = D
            Z = g4.integrate(sd, "none", 0.0, xx, v, np.random.default_rng(31 + 7 * sd + SPLIT_ID[s]))
            Fl[s], Fs[s] = run(sd, Z, xx[-1] + D, np.random.default_rng(6000 + 13 * sd + SPLIT_ID[s]))
        per[sd] = {"avg_long": g7.best_predict(Fl, data)[1], "avg_short": g7.best_predict(Fs, data)[1]}
        print(sd, per[sd], flush=True)
    S = {k: [float(np.mean([per[s][k] for s in SEEDS])), float(np.std([per[s][k] for s in SEEDS]))] for k in per[SEEDS[0]]}
    print(S)
    os.makedirs("results", exist_ok=True)
    json.dump({"amp": AMP, "summary_mean_sd": S, "per_seed": {str(k): v for k, v in per.items()}},
              open("results/gate8b_receipt.json", "w"), indent=2)


if __name__ == "__main__":
    main()
