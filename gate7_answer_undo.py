"""Gate 7 - can the ANSWER itself tell you how to undo the read?

Sol's correction to Gate 6b: exact restoration need not need a saved pre-query copy; known,
invertible dynamics can sometimes reconstruct it. The cheapest candidate: the decoded answer.
To undo a ping at relative phase alpha_k you kick at relative phase -alpha_k. If alpha_k could be
computed from the decoded goal vector D_hat (two numbers) and the bank's known structure
(b_k, d_k), the undo needs no per-unit copy at all: alpha_hat_k = b_k d_k . D_hat.

Fairness note: the Gate 6 one-bank listener reads every unit's phase BEFORE pinging, so it
already holds a per-unit pre-query copy. This gate uses that copy only as a labelled reference;
the question is whether the 2-number answer can replace it.

Protocol (Gates 2/4/6 bank: 40 units, beta = 0, sigma = 0.05; goal |D| <= 0.35):
  window 1: ping +a e^{i phi}, listen 16 substeps; a reader (ridge / kNN, chosen on validation)
            decodes D_hat from window-1 one-bank features (as Gate 6). The answer is identical for
            every undo scheme, so no accuracy matching is needed.
  undo (test episodes), all kicks of size a*sqrt(mu), u = the unit's current phase:
    none            -
    sign_flip       -a e^{i phi}                    (memoryless; Gate 6)
    answer          a u e^{-i b d.D_hat}            (built from the 2-number decoded answer)
    answer_true     a u e^{-i b d.D}                (same, with the TRUE D: isolates answer error
                                                     from what any 2-number answer cannot carry)
    prequery        a u e^{-i(phi - theta_0)}       (per-unit pre-query copy; exact, Gate 6b)
  then 32 substeps of relaxation; damage = RMS per-unit phase difference against the unqueried twin
  (co-rotating shared noise). a in {0.1, 0.2, 0.3}.

Kill conditions, fixed before running:
  K7a  `answer` leaves <= 0.5x the sign_flip damage at a = 0.2. Otherwise the 2-number answer
       cannot replace a per-unit copy.
  K7b  diagnosis: if K7a fails, `answer_true` also fails it (the limit is information the answer
       does not carry, not decoding error).
Prediction written with the kill conditions: K7a fails. Each unit's stored phase carries its own
noise drift (0.71 rad RMS, Gate 4), which a 2-number position estimate averages away, so
alpha_hat is wrong by about that much per unit and the residual is first order in a.
"""
import json
import os
import time

import numpy as np

import gate2_path as g2
import gate4_ping as g4
import gate5_sweeps as g5
import gate6_restore as g6

MU = g2.MU
E, OFF, SEEDS = g6.E, g6.OFF, g6.SEEDS
AMPS = (0.1, 0.2, 0.3)
SCHEMES = ("none", "sign_flip", "answer", "answer_true", "prequery")
SPLIT_ID = {"tr": 1, "va": 2, "te": 3}


def best_predict(F, Y):
    """Return (val-selected) test predictions and the test error."""
    best = None
    for f in (g4.ridge, g4.knn):
        pv = f(F["tr"], Y["tr"], F["va"], Y["va"], F["va"])
        ev = np.linalg.norm(pv - Y["va"], axis=1).mean()
        if best is None or ev < best[0]:
            best = (ev, f)
    P = best[1](F["tr"], Y["tr"], F["va"], Y["va"], F["te"])
    return P, float(np.linalg.norm(P - Y["te"], axis=1).mean())


def window1(sd, Z, q, amp, rng):
    _, d, b = g2.bank(sd)
    phi = b[:, None] * (d @ q.T)
    zp, zu = Z.copy(), Z.copy()
    th0 = np.angle(zp)
    zp = zp + amp * np.sqrt(MU) * np.exp(1j * phi)
    feats = []
    for s in range(1, g6.H + 1):
        zp, zu = g6.step(zp, zu, rng)
        if s in g6.SAMPLE:
            feats.append(np.abs(zp) - np.sqrt(MU)); feats.append(np.angle(zp * np.exp(-1j * th0)))
    return np.concatenate(feats, 0).T, zp, zu, phi, th0


def undo_damage(sd, zp, zu, phi, th0, D_use, amp, scheme, seed):
    _, d, b = g2.bank(sd)
    u = zp / np.abs(zp)
    k = amp * np.sqrt(MU)
    if scheme == "sign_flip":
        zp = zp - k * np.exp(1j * phi)
    elif scheme in ("answer", "answer_true"):
        zp = zp + k * u * np.exp(-1j * b[:, None] * (d @ D_use.T))
    elif scheme == "prequery":
        zp = zp + k * u * np.exp(-1j * (phi - th0))
    rng = np.random.default_rng(seed)
    for _ in range(2 * g6.RELAX):
        zp, zu = g6.step(zp, zu, rng)
    return float(np.sqrt((np.angle(zp * np.conj(zu)) ** 2).mean()))


def run_seed(sd):
    data, Z = {}, {}
    for s in E:
        xx, v = g2.trajectories(np.random.default_rng(OFF[s] + sd), E[s])
        rng = np.random.default_rng(900 + OFF[s] + sd)
        ang = rng.uniform(0, 2 * np.pi, E[s]); rad = 0.35 * np.sqrt(rng.uniform(0, 1, E[s]))
        D = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
        data[s] = {"D": D, "q": xx[-1] + D}
        Z[s] = g4.integrate(sd, "none", 0.0, xx, v, np.random.default_rng(31 + 7 * sd + SPLIT_ID[s]))
    Y = {s: data[s]["D"] for s in E}
    res = {}
    for amp in AMPS:
        F, st = {}, {}
        for s in E:
            out = window1(sd, Z[s], data[s]["q"], amp, np.random.default_rng(5000 + 13 * sd + SPLIT_ID[s]))
            F[s], st[s] = out[0], out[1:]
        P, err = best_predict(F, Y)
        zp, zu, phi, th0 = st["te"]
        row = {"answer_error": err}
        for sch in SCHEMES:
            D_use = P if sch == "answer" else Y["te"]
            row[sch] = undo_damage(sd, zp, zu, phi, th0, D_use, amp, sch, 7000 + sd)
        res[f"a{amp:g}"] = row
        print(sd, amp, {k: round(v, 4) for k, v in row.items()}, flush=True)
    return res


def main():
    t0 = time.time()
    per = {sd: run_seed(sd) for sd in SEEDS}
    S = {a: {m: [float(np.mean([per[s][a][m] for s in SEEDS])), float(np.std([per[s][a][m] for s in SEEDS]))]
             for m in per[SEEDS[0]][a]} for a in per[SEEDS[0]]}
    for a, row in S.items():
        print(a, {m: round(v[0], 4) for m, v in row.items()})
    r = S["a0.2"]
    V = {"K7a_answer_halves_signflip": bool(r["answer"][0] <= 0.5 * r["sign_flip"][0]),
         "K7a_values": {"answer": r["answer"][0], "sign_flip": r["sign_flip"][0]},
         "K7b_true_answer_also_fails": bool(r["answer_true"][0] > 0.5 * r["sign_flip"][0]),
         "K7b_values": {"answer_true": r["answer_true"][0], "prequery": r["prequery"][0], "none": r["none"][0]}}
    out = {"setup": {"amps": AMPS, "schemes": SCHEMES, "episodes": E, "seeds": SEEDS,
                     "listen_substeps": g6.H, "relax_substeps": 2 * g6.RELAX},
           "summary_mean_sd": S, "verdicts": V, "per_seed": {str(k): v for k, v in per.items()},
           "seconds": round(time.time() - t0, 1)}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/gate7_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), out["seconds"])


if __name__ == "__main__":
    main()
