"""Gate 4 - ping the memory with a goal, listen for the vector to it.

Prompted by Sol's correction to Gate 3: equal singular values do not mean equal responses.
An uncoupled oscillator's response operator ROTATES with its stored phase, so a fixed probe
gets a phase-dependent answer with no coupling at all. Gate 3 only measured rotation-invariant
change. This gate measures actual probe answers on a task whose answer depends jointly on the
stored state and the probe, along the real noisy trajectory (Sol's two qualifications).

Task (goal-vector readout):
  1. Path integration exactly as Gate 2 (40 units, beta = 0, w0 = 0, noise sigma = 0.05,
     200 steps, phases set at a start landmark, then velocity only).
  2. At the end, a goal q near the true position x (q = x + D, |D| <= 0.35) is given as a
     PING: an instantaneous kick a*sqrt(mu)*exp(i phi_k), with phi_k = b_k d_k . q, the phase
     pattern the bank would hold at q.
  3. LISTEN: evolve H = 4 time units standing still (v = 0), with a probed and an unprobed copy
     sharing the same noise. Features: each unit's phase advance arg(z_probed conj(z_unprobed))
     at t = 1, 2, 4. For an uncoupled unit, to first order, it is (a) sin(phi_k - theta_k):
     interference between the goal and the stored phase.
  4. A reader (ridge, or kNN; hyperparameters on validation) maps the features to D = q - x.

Conditions (coupling during integration / during the listen window):
  none         0 / 0
  fixed_vortex 0.03 / 0.03     (steady coupling, costs integration accuracy, Gate 3)
  gated_vortex 0 / k_on        k_on in {0.1, 0.3}: coupling switched on by the ping only
  gated_linear 0 / k_on        same with ordinary coupling
  Ping amplitude a in {0.05, 0.3}.
  ORACLE: the same readers on sin/cos(phi_k - theta_k) computed from the actual noisy end
  phases. This is what a perfect, noiseless, non-disturbing ping would deliver.

Also measured: READ DAMAGE. A ridge reader trained on end phases -> position is applied to the
probed copy after the listen window; damage = its error minus the error before the ping.

Kill conditions, fixed before running:
  K1  coupling unnecessary: uncoupled ping, best reader, D error <= 0.5 x the error of always
      answering D = 0. Otherwise the stored phase does not usefully change probe answers.
  K2  physical ping is near ideal: the best uncoupled ping condition is within 1.25x the
      ORACLE error. Otherwise pinging loses most of what is stored.
  K3  vortex earns its place: some vortex condition beats the best uncoupled ping by >= 10%
      on D error, with the gap larger than 2 seed-sd. Otherwise vortex adds nothing here.
"""
import json
import os
import time

import numpy as np

import gate2_path as g2

N, MU, SIGMA, DT, SUB = g2.N, g2.MU, g2.SIGMA, g2.DT, g2.SUB
E_TR, E_VA, E_TE = 1000, 300, 500
H_SUB = 16                       # listen window: 16 substeps = 4 time units
SAMPLE = (4, 8, 16)              # substeps at which phase advance is read (t = 1, 2, 4)
RADIUS = 0.35
SEEDS = (10, 11, 12)
AMPS = (0.05, 0.3)
CONDS = [("none", "none", 0.0, 0.0), ("fixed_vortex", "vortex", 0.03, 0.03),
         ("gated_vortex_0.1", "vortex", 0.0, 0.1), ("gated_vortex_0.3", "vortex", 0.0, 0.3),
         ("gated_linear_0.1", "linear", 0.0, 0.1), ("gated_linear_0.3", "linear", 0.0, 0.3)]


def rhs(z, lin, M, mode, kappa):
    out = lin * z - (z.real ** 2 + z.imag ** 2) * z
    if kappa > 0:
        out = out + kappa * (M @ np.conj(z) if mode == "vortex" else M @ z)
    return out


def rk4(z, lin, M, mode, kappa):
    k1 = rhs(z, lin, M, mode, kappa); k2 = rhs(z + 0.5 * DT * k1, lin, M, mode, kappa)
    k3 = rhs(z + 0.5 * DT * k2, lin, M, mode, kappa); k4 = rhs(z + DT * k3, lin, M, mode, kappa)
    return z + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def noise(rng, shape):
    return SIGMA * np.sqrt(DT) * (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / np.sqrt(2)


def integrate(seed, mode, kappa, x, v, rng):
    M, d, b = g2.bank(seed)
    z = np.sqrt(MU) * np.exp(1j * b[:, None] * (d @ x[0].T))
    for t in range(g2.STEPS):
        lin = MU + 1j * (b[:, None] * (d @ v[t].T) / (SUB * DT))
        for _ in range(SUB):
            z = rk4(z, lin, M, mode, kappa) + noise(rng, z.shape)
    return z


def listen(seed, z_end, q, amp, mode, kappa, rng):
    M, d, b = g2.bank(seed)
    phi = b[:, None] * (d @ q.T)
    zp = z_end + amp * np.sqrt(MU) * np.exp(1j * phi)
    zu = z_end.copy()
    feats = []
    for s in range(1, H_SUB + 1):
        n = noise(rng, z_end.shape)
        zp = rk4(zp, MU, M, mode, kappa) + n
        zu = rk4(zu, MU, M, mode, kappa) + n
        if s in SAMPLE:
            feats.append(np.angle(zp * np.conj(zu)))
    return np.concatenate(feats, 0).T, zp


def phase_feats(z):
    u = z / np.maximum(np.abs(z), 1e-9)
    return np.concatenate([u.real, u.imag], 0).T


def oracle_feats(seed, z_end, q):
    M, d, b = g2.bank(seed)
    rel = b[:, None] * (d @ q.T) - np.angle(z_end)
    return np.concatenate([np.sin(rel), np.cos(rel)], 0).T


# ---------- readers -------------------------------------------------------------------------
def standardize(Ftr, *others):
    m, s = Ftr.mean(0), Ftr.std(0) + 1e-9
    return [(F - m) / s for F in (Ftr,) + others]


def ridge(Ftr, Ytr, Fva, Yva, Fte):
    Ftr, Fva, Fte = standardize(Ftr, Fva, Fte)
    add = lambda F: np.concatenate([F, np.ones((len(F), 1))], 1)
    best = None
    for lam in (1e-2, 1e-1, 1.0, 10.0, 100.0):
        W = np.linalg.solve(add(Ftr).T @ add(Ftr) + lam * np.eye(Ftr.shape[1] + 1), add(Ftr).T @ Ytr)
        e = np.linalg.norm(add(Fva) @ W - Yva, axis=1).mean()
        if best is None or e < best[0]:
            best = (e, W)
    return add(Fte) @ best[1]


def knn(Ftr, Ytr, Fva, Yva, Fte):
    Ftr, Fva, Fte = standardize(Ftr, Fva, Fte)

    def pred(F, k):
        d2 = (F ** 2).sum(1)[:, None] - 2 * F @ Ftr.T + (Ftr ** 2).sum(1)[None]
        idx = np.argpartition(d2, k, axis=1)[:, :k]
        return Ytr[idx].mean(1)
    k = min((5, 15, 40), key=lambda k: np.linalg.norm(pred(Fva, k) - Yva, axis=1).mean())
    return pred(Fte, k)


def score(F, Y, readers=("ridge", "knn")):
    out = {}
    for r in readers:
        P = (ridge if r == "ridge" else knn)(F["tr"], Y["tr"], F["va"], Y["va"], F["te"])
        out[r] = float(np.linalg.norm(P - Y["te"], axis=1).mean())
    return out


# ---------- experiment ----------------------------------------------------------------------
def run_seed(sd):
    split = {"tr": (100, E_TR), "va": (200, E_VA), "te": (300, E_TE)}
    data = {}
    for k, (off, E) in split.items():
        x, v = g2.trajectories(np.random.default_rng(off + sd), E)
        rng = np.random.default_rng(900 + off + sd)
        ang = rng.uniform(0, 2 * np.pi, E); rad = RADIUS * np.sqrt(rng.uniform(0, 1, E))
        D = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
        data[k] = {"x": x[-1], "D": D, "q": x[-1] + D, "x_all": x, "v": v}

    res = {"predict_zero_error": float(np.linalg.norm(data["te"]["D"], axis=1).mean())}
    ends = {}
    for name, mode, k_int, _ in CONDS:
        key = (mode if k_int > 0 else "none", k_int)
        if key not in ends:
            ends[key] = {s: integrate(sd, key[0], k_int, data[s]["x_all"], data[s]["v"],
                                      np.random.default_rng(31 + 7 * sd + {"tr": 1, "va": 2, "te": 3}[s])) for s in split}
    Y = {s: data[s]["D"] for s in split}
    X = {s: data[s]["x"] for s in split}

    for key, Z in ends.items():
        Fo = {s: oracle_feats(sd, Z[s], data[s]["q"]) for s in split}
        res[f"oracle_int_{key[0]}_{key[1]:g}"] = score(Fo, Y)
        Fp = {s: phase_feats(Z[s]) for s in split}
        # absolute position straight from the phases (added after the first run, for a fair comparison)
        res[f"position_before_ping_int_{key[0]}_{key[1]:g}"] = score(Fp, X)

    for name, mode, k_int, k_read in CONDS:
        Z = ends[(mode if k_int > 0 else "none", k_int)]
        Fp_tr = phase_feats(Z["tr"])
        for amp in AMPS:
            F, post = {}, {}
            for s in split:
                F[s], post[s] = listen(sd, Z[s], data[s]["q"], amp, mode if k_read > 0 else "none", k_read,
                                       np.random.default_rng(77 + 11 * sd + {"tr": 1, "va": 2, "te": 3}[s]))
            sc = score(F, Y)
            # read damage: position reader trained on pre-ping phases, applied after the window
            Fpos = {"tr": Fp_tr, "va": phase_feats(Z["va"]), "te": phase_feats(post["te"])}
            after = score(Fpos, X, ("knn",))["knn"]
            before = res[f"position_before_ping_int_{mode if k_int > 0 else 'none'}_{k_int:g}"]["knn"]
            res[f"{name}_a{amp:g}"] = {**sc, "best": min(sc.values()),
                                       "position_error_before_knn": before, "position_error_after_knn": after}
    return res


def main():
    t0 = time.time()
    per = {}
    for sd in SEEDS:
        per[sd] = run_seed(sd)
        print("seed", sd, f"{time.time() - t0:.0f}s", flush=True)
    keys = [k for k in per[SEEDS[0]] if isinstance(per[SEEDS[0]][k], dict)]
    summary = {"predict_zero_error": float(np.mean([per[s]["predict_zero_error"] for s in SEEDS]))}
    for k in keys:
        summary[k] = {m: [float(np.mean([per[s][k][m] for s in SEEDS])), float(np.std([per[s][k][m] for s in SEEDS]))]
                      for m in per[SEEDS[0]][k]}
    for k, v in summary.items():
        print(k, v if not isinstance(v, dict) else {m: round(x[0], 4) for m, x in v.items()})

    S = summary
    zero = S["predict_zero_error"]
    unc = {k: S[k]["best"] for k in S if k.startswith("none_a")}
    best_unc_key = min(unc, key=unc.get)
    best_unc, best_unc_sd = S[best_unc_key]["best"]
    oracle = min(S["oracle_int_none_0"].values(), key=lambda q: q[0])[0]
    vort = {k: S[k]["best"] for k in S if "vortex" in k and isinstance(S[k], dict) and "best" in S[k]}
    best_v_key = min(vort, key=lambda k: vort[k][0])
    bv, bv_sd = vort[best_v_key]
    V = {"K1_uncoupled_ping_beats_half_of_zero": bool(best_unc <= 0.5 * zero),
         "K1_values": [best_unc, zero],
         "K2_uncoupled_within_1.25x_oracle": bool(best_unc <= 1.25 * oracle),
         "K2_values": [best_unc, oracle],
         "K3_vortex_beats_uncoupled_by_10pct": bool(bv <= 0.9 * best_unc and best_unc - bv > 2 * max(bv_sd, best_unc_sd)),
         "K3_values": {"best_vortex": [best_v_key, bv], "best_uncoupled": [best_unc_key, best_unc]}}
    out = {"setup": {"N": N, "episodes": [E_TR, E_VA, E_TE], "listen_substeps": H_SUB, "sample_substeps": SAMPLE,
                     "goal_radius": RADIUS, "amps": AMPS, "conditions": CONDS, "seeds": SEEDS},
           "summary_mean_sd": summary, "per_seed": {str(k): v for k, v in per.items()}, "verdicts": V,
           "seconds": round(time.time() - t0, 1)}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/gate4_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), out["seconds"])


if __name__ == "__main__":
    main()
