"""Gate 8 - can the bank be queried and restored through a cheap interface?

Sol's open hurdle: VMN and VMNClaude read every unit through a large decoder. A physical memory
read by pinging is interesting only if a few wires suffice. Two cheap listeners are compared
with the full per-unit listener (Gates 6/7), all on the Gates 2/4 bank (40 units, beta = 0,
sigma = 0.05; goal |D| <= 0.35; ping strength a):

  full       every unit's amplitude and phase change, 16 substeps (4 time units) - reference
  chan_K     K in {1, 2, 4} summed channels y_c(t) = sum_k w_ck z_k(t) / sqrt(N), fixed random
             unit-modulus complex weights, all units at the same frequency; features: y_c just
             before the ping and its change at t = 0.25, 0.5, 1, 2, 4 (Sol's aggregate listener,
             in this bank)
  wire_fdm   ONE wire. Each unit rests at its own frequency, omega_k = (k - N/2) * dw, with
             dw = 2 pi / (L dt), L = 48 substeps (12 time units), so the frequencies are exactly
             orthogonal over a window. y(t) = sum_k z_k(t) e^{i omega_k t} + measurement noise.
             The wire is listened to for one window before the ping and one after; demodulating
             each window at omega_k gives every unit's averaged state, and the features are each
             unit's phase and amplitude change between the two windows. Measurement noise
             sigma_m in {0, 0.5, 2} per complex sample (the wire's own amplitude is ~ sqrt(N mu)
             = 4.5). This is frequency-multiplexed readout, the known way to read many
             resonators on one line (NEEDS.md); the gate asks whether it carries this task.
             Isochronous units: a resting frequency is a rotation of the frame and does not
             change the dynamics, so the bank is simulated in each unit's rotating frame.

For wire_fdm the read is followed by the sign-flip undo (Gate 6) and the damage is measured
against an unqueried twin after relaxation; `none` is the same protocol without the undo.
Readers: ridge / kNN, best on validation. a in {0.2, 0.3}; 3 seeds; 1000 / 300 / 500 episodes.
Injection is idealised in every condition (each unit receives its kick); with frequency
multiplexing the same pattern could in principle be sent down the one wire as a sum of tones.

Kill conditions, fixed before running (at a = 0.2):
  K8a  some chan_K reaches answer error <= 1.5x the full listener's.
  K8b  wire_fdm with sigma_m = 0.5 reaches answer error <= 1.5x the full listener's.
  K8c  within the wire_fdm protocol, the sign-flip undo leaves <= 1/3 of the no-undo damage.
Prediction: K8a fails (Sol's four-channel result), K8b passes (frequency multiplexing is a
per-unit readout in disguise, paid for with 6x the listening time), K8c passes (the undo needs
no readout at all).
"""
import json
import os
import time

import numpy as np

import gate2_path as g2
import gate4_ping as g4
import gate6_restore as g6
import gate7_answer_undo as g7

N, MU, DT = g2.N, g2.MU, g2.DT
E, OFF, SEEDS, SPLIT_ID = g6.E, g6.OFF, g6.SEEDS, g7.SPLIT_ID
AMPS = (0.2, 0.3)
KS = (1, 2, 4)
L = 48
DW = 2 * np.pi / (L * DT)
OMEGA = (np.arange(N) - N / 2) * DW
SIG_M = (0.0, 0.5, 2.0)


def score(F, Y):
    return g7.best_predict(F, Y)[1]


def short_window(sd, Z, q, amp, rng, W):
    """Gate 6/7 window: full features and K-channel features from the same run."""
    _, d, b = g2.bank(sd)
    phi = b[:, None] * (d @ q.T)
    zp, zu = Z.copy(), Z.copy()
    th0 = np.angle(zp)
    y0 = {K: (W[K] @ zp) / np.sqrt(N) for K in KS}
    zp = zp + amp * np.sqrt(MU) * np.exp(1j * phi)
    full, ch = [], {K: [y0[K].real, y0[K].imag] for K in KS}
    for s in range(1, g6.H + 1):
        zp, zu = g6.step(zp, zu, rng)
        if s in g6.SAMPLE:
            full.append(np.abs(zp) - np.sqrt(MU)); full.append(np.angle(zp * np.exp(-1j * th0)))
            for K in KS:
                dy = (W[K] @ zp) / np.sqrt(N) - y0[K]
                ch[K] += [dy.real, dy.imag]
    return np.concatenate(full, 0).T, {K: np.concatenate(ch[K], 0).T for K in KS}


def fdm(sd, Z, q, amp, rng, mrng, undo):
    """One wire, frequency-multiplexed. Returns features per noise level and the damage."""
    _, d, b = g2.bank(sd)
    phi = b[:, None] * (d @ q.T)
    zp, zu = Z.copy(), Z.copy()
    Ecount = Z.shape[1]
    noise = {s: mrng.normal(size=(2 * L, Ecount)) + 1j * mrng.normal(size=(2 * L, Ecount)) for s in SIG_M}
    demod = {s: [np.zeros((N, Ecount), complex), np.zeros((N, Ecount), complex)] for s in SIG_M}

    def record(n, w):
        rot = np.exp(1j * OMEGA * n * DT)[:, None]
        y = (zp * rot).sum(0)
        for s in SIG_M:
            demod[s][w] += (y + s / np.sqrt(2) * noise[s][n]) [None, :] * np.conj(rot)

    for n in range(L):
        record(n, 0)
        zp, zu = g6.step(zp, zu, rng)
    zp = zp + amp * np.sqrt(MU) * np.exp(1j * phi)
    for n in range(L, 2 * L):
        zp, zu = g6.step(zp, zu, rng)
        record(n, 1)
    feats = {}
    for s in SIG_M:
        pre, post = demod[s][0] / L, demod[s][1] / L
        feats[s] = np.concatenate([np.angle(post * np.conj(pre)), np.abs(post) - np.abs(pre)], 0).T
    if undo:
        zp = zp - amp * np.sqrt(MU) * np.exp(1j * phi)
    for _ in range(2 * g6.RELAX):
        zp, zu = g6.step(zp, zu, rng)
    return feats, float(np.sqrt((np.angle(zp * np.conj(zu)) ** 2).mean()))


def run_seed(sd):
    wr = np.random.default_rng(4000 + sd)
    W = {K: np.exp(1j * wr.uniform(0, 2 * np.pi, (K, N))) for K in KS}
    data, Z = {}, {}
    for s in E:
        xx, v = g2.trajectories(np.random.default_rng(OFF[s] + sd), E[s])
        rng = np.random.default_rng(900 + OFF[s] + sd)
        ang = rng.uniform(0, 2 * np.pi, E[s]); rad = 0.35 * np.sqrt(rng.uniform(0, 1, E[s]))
        D = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
        data[s] = {"D": D, "q": xx[-1] + D}
        Z[s] = g4.integrate(sd, "none", 0.0, xx, v, np.random.default_rng(31 + 7 * sd + SPLIT_ID[s]))
    Y = {s: data[s]["D"] for s in E}
    res = {"predict_zero_error": float(np.linalg.norm(Y["te"], axis=1).mean())}
    for amp in AMPS:
        Ff, Fc, Fw, dmg = {}, {K: {} for K in KS}, {s: {} for s in SIG_M}, {}
        for sp in E:
            f, c = short_window(sd, Z[sp], data[sp]["q"], amp, np.random.default_rng(5000 + 13 * sd + SPLIT_ID[sp]), W)
            Ff[sp] = f
            for K in KS:
                Fc[K][sp] = c[K]
            for undo in ((False, True) if sp == "te" else (True,)):
                fw, dm = fdm(sd, Z[sp], data[sp]["q"], amp, np.random.default_rng(6000 + 13 * sd + SPLIT_ID[sp]),
                             np.random.default_rng(8000 + 13 * sd + SPLIT_ID[sp]), undo)
                if sp == "te":
                    dmg["fdm_undo" if undo else "fdm_none"] = dm
            for s in SIG_M:
                Fw[s][sp] = fw[s]
        row = {"full": score(Ff, Y)}
        for K in KS:
            row[f"chan_{K}"] = score(Fc[K], Y)
        for s in SIG_M:
            row[f"wire_fdm_sigma{s:g}"] = score(Fw[s], Y)
        row.update({f"damage_{k}": v for k, v in dmg.items()})
        res[f"a{amp:g}"] = row
        print(sd, amp, {k: round(v, 4) for k, v in row.items()}, flush=True)
    return res


def main():
    t0 = time.time()
    per = {sd: run_seed(sd) for sd in SEEDS}
    S = {"predict_zero_error": float(np.mean([per[s]["predict_zero_error"] for s in SEEDS]))}
    for a in [k for k in per[SEEDS[0]] if k.startswith("a")]:
        S[a] = {m: [float(np.mean([per[s][a][m] for s in SEEDS])), float(np.std([per[s][a][m] for s in SEEDS]))]
                for m in per[SEEDS[0]][a]}
        print(a, {m: round(v[0], 4) for m, v in S[a].items()})
    r = S["a0.2"]
    full = r["full"][0]
    best_k = min(KS, key=lambda K: r[f"chan_{K}"][0])
    V = {"K8a_few_channels_within_1.5x": bool(r[f"chan_{best_k}"][0] <= 1.5 * full),
         "K8a_values": {"full": full, **{f"chan_{K}": r[f"chan_{K}"][0] for K in KS}},
         "K8b_one_fdm_wire_within_1.5x": bool(r["wire_fdm_sigma0.5"][0] <= 1.5 * full),
         "K8b_values": {s: r[f"wire_fdm_sigma{s:g}"][0] for s in SIG_M},
         "K8c_undo_within_fdm_cuts_3x": bool(r["damage_fdm_undo"][0] <= r["damage_fdm_none"][0] / 3),
         "K8c_values": {"none": r["damage_fdm_none"][0], "undo": r["damage_fdm_undo"][0]}}
    out = {"setup": {"amps": AMPS, "channels": KS, "fdm_window_substeps": L, "fdm_dw": DW, "sigma_m": SIG_M,
                     "episodes": E, "seeds": SEEDS},
           "summary_mean_sd": S, "verdicts": V, "per_seed": {str(k): v for k, v in per.items()},
           "seconds": round(time.time() - t0, 1)}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/gate8_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), out["seconds"])


if __name__ == "__main__":
    main()
