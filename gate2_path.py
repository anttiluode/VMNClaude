"""Gate 2 - path integration with velocity-controlled oscillators (entorhinal setting).

Why this task: in the oscillatory-interference model of grid cells, each oscillator's
frequency relative to a baseline is b_k (d_k . v): it turns one way when the animal
moves along d_k and the other way when it moves against it. Both rotation senses are
built in, which is where Gate 1b found the vortex (conjugate) coupling switches on.

Units (rotating frame of a baseline rhythm w0, read out by interference with it):

    dz_k/dt = (mu + i(w0 + b_k d_k . v)) z_k - (1 + i beta)|z_k|^2 z_k + kappa C(z)_k + sigma dW_k

    C(z) = M conj(z)  vortex      C(z) = M z  linear      C(z) = 0  none

M is the Gate 1 vortex matrix (fixed random unit positions, spectral norm 1).
Each episode starts with phases set to the true start position (a landmark at t = 0),
then only velocity comes in; independent noise sigma on every unit makes phases drift.
Readout: one ridge regression from phase features [cos, sin of arg(z_k e^{-i w0 t})]
to position, trained on all time steps of training episodes.

Baseline: echo state network (2N tanh units, leaky) given velocity plus a start-position
pulse at t = 0, with the same noise per state variable.

Kill conditions, fixed before running:
  K1  with w0 = 0 (both senses present), vortex coupling must beat BOTH linear coupling
      and no coupling on test end-of-episode position error; otherwise the conjugate
      coupling does not help path integration.
  K2  the best VMN variant must beat the ESN; otherwise not useful here.
  K3  mechanism control: with w0 = 3 (every unit turns the same way in the lab frame),
      the vortex coupling's advantage over no coupling must shrink by at least half;
      otherwise the rotation-sense explanation is wrong.
"""
import json
import os
import time

import numpy as np

import gate1_layer as g1

N = g1.N                      # 40 units
STEPS, SUB, DT = 200, 4, 0.25
E_TR, E_VA, E_TE = 300, 100, 200
MU = 0.5
SIGMA = 0.05
SCALES = (2 * np.pi * 1.0, 2 * np.pi * 2.0, 2 * np.pi * 4.0)   # radians per unit distance


def trajectories(rng, E):
    """Smooth random walks in the unit box with reflecting walls. Returns x (T+1,E,2), v (T,E,2)."""
    x = np.zeros((STEPS + 1, E, 2))
    v = np.zeros((STEPS, E, 2))
    x[0] = rng.uniform(0.1, 0.9, (E, 2))
    vel = rng.normal(size=(E, 2)) * 0.02
    for t in range(STEPS):
        vel = 0.9 * vel + 0.1 * rng.normal(size=(E, 2)) * 0.03
        nx = x[t] + vel
        for d in range(2):
            lo, hi = nx[:, d] < 0, nx[:, d] > 1
            nx[lo, d] = -nx[lo, d]
            nx[hi, d] = 2 - nx[hi, d]
            vel[lo | hi, d] *= -1
        v[t] = nx - x[t]
        x[t + 1] = nx
    return x, v


def bank(seed):
    rng = np.random.default_rng(1000 + seed)
    M, _ = g1.vortex_matrix(rng)
    ang = rng.uniform(0, 2 * np.pi, N)
    d = np.stack([np.cos(ang), np.sin(ang)], 1)
    b = np.array([SCALES[k % 3] for k in range(N)])
    return M, d, b


def run_vmn(seed, mode, kappa, w0, beta, x, v, noise_seed, sigma=SIGMA):
    M, d, b = bank(seed)
    E = x.shape[1]
    rng = np.random.default_rng(noise_seed)
    r0 = np.sqrt(MU)
    z = r0 * np.exp(1j * b[:, None] * (d @ x[0].T))          # phases set to true start
    feats = np.empty((STEPS + 1, E, 2 * N))
    t_now = 0.0

    def feat(z, t):
        u = z * np.exp(-1j * (w0 - beta * MU) * t)   # resting cycle turns at w0 - beta*mu
        ph = u / np.maximum(np.abs(u), 1e-9)
        return np.concatenate([ph.real, ph.imag], 0).T

    feats[0] = feat(z, 0.0)
    for t in range(STEPS):
        freq = w0 + b[:, None] * (d @ v[t].T) / (SUB * DT)   # phase advance b d.dx over the step
        lin = MU + 1j * freq

        def f(z):
            if mode == "vortex":
                c = M @ np.conj(z)
            elif mode == "linear":
                c = M @ z
            else:
                c = 0.0
            return lin * z - (1 + 1j * beta) * (z.real ** 2 + z.imag ** 2) * z + kappa * c

        for _ in range(SUB):
            k1 = f(z); k2 = f(z + 0.5 * DT * k1)
            k3 = f(z + 0.5 * DT * k2); k4 = f(z + DT * k3)
            z = z + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            z = z + sigma * np.sqrt(DT) * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape)) / np.sqrt(2)
            t_now += DT
        feats[t + 1] = feat(z, t_now)
    return feats


def run_esn(seed, rho, s, leak, x, v, noise_seed, sigma=SIGMA):
    rng = np.random.default_rng(2000 + seed)
    n = 2 * N
    W = rng.normal(size=(n, n))
    W *= rho / np.max(np.abs(np.linalg.eigvals(W)))
    Win = rng.uniform(-1, 1, (n, 4))
    nrng = np.random.default_rng(noise_seed)
    E = x.shape[1]
    h = np.zeros((n, E))
    feats = np.empty((STEPS + 1, E, n))
    # start pulse carries the true start position
    u0 = np.concatenate([np.zeros((2, E)), x[0].T - 0.5], 0)
    h = np.tanh(W @ h + s * Win @ u0)
    feats[0] = h.T
    for t in range(STEPS):
        u = np.concatenate([v[t].T / 0.02, np.zeros((2, E))], 0)
        h = (1 - leak) * h + leak * np.tanh(W @ h + s * Win @ u)
        h = h + sigma * nrng.normal(size=h.shape)
        feats[t + 1] = h.T
    return feats


def fit_and_score(Ftr, xtr, Fev, xev):
    """Ridge from features to position over all steps; lambda chosen on the eval split it is given."""
    X = Ftr.reshape(-1, Ftr.shape[-1]); X = np.concatenate([X, np.ones((len(X), 1))], 1)
    Y = xtr.reshape(-1, 2)
    Xe = Fev.reshape(-1, Fev.shape[-1]); Xe = np.concatenate([Xe, np.ones((len(Xe), 1))], 1)
    best = None
    for lam in (1e-4, 1e-2, 1.0, 1e2):
        Wt = np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ Y)
        P = (Xe @ Wt).reshape(xev.shape)
        err = np.linalg.norm(P - xev, axis=-1)                # (T+1, E)
        score = err[-1].mean()
        if best is None or score < best[0]:
            best = (score, err.mean(), lam, Wt)
    return best


def evaluate(make_feats, seeds, tune_only=False):
    """Train readout on train episodes, pick lambda on val; report test (unless tune_only)."""
    out = []
    for sd in seeds:
        xtr, vtr = trajectories(np.random.default_rng(100 + sd), E_TR)
        xva, vva = trajectories(np.random.default_rng(200 + sd), E_VA)
        xte, vte = trajectories(np.random.default_rng(300 + sd), E_TE)
        Ftr = make_feats(sd, xtr, vtr, 11 + 7 * sd)
        Fva = make_feats(sd, xva, vva, 12 + 7 * sd)
        val_end, _, lam, Wt = fit_and_score(Ftr, xtr, Fva, xva)
        if tune_only:
            out.append({"val_end_error": float(val_end)})
            continue
        Fte = make_feats(sd, xte, vte, 13 + 7 * sd)
        Xe = Fte.reshape(-1, Fte.shape[-1]); Xe = np.concatenate([Xe, np.ones((len(Xe), 1))], 1)
        P = (Xe @ Wt).reshape(xte.shape)
        err = np.linalg.norm(P - xte, axis=-1)
        out.append({"val_end_error": float(val_end), "test_end_error": float(err[-1].mean()),
                    "test_mean_error": float(err.mean()),
                    "test_error_curve": [float(e) for e in err.mean(1)[::20]]})
    return out


def main():
    t0 = time.time()
    TUNE, TEST = (0, 1), (10, 11, 12)
    res = {"setup": {"N": N, "steps": STEPS, "substeps": SUB, "dt": DT, "mu": MU, "sigma": SIGMA,
                     "scales_rad_per_unit": SCALES, "episodes": [E_TR, E_VA, E_TE],
                     "error_if_always_predicting_box_centre": 0.3826},
           "vmn": {}, "esn": {}}
    for w0 in (0.0, 3.0):
        for beta in (0.0, 2.0):
            for mode in ("vortex", "linear", "none"):
                key = f"w0{w0:g}_{mode}_beta{beta:g}"
                grid = []
                for kappa in ((0.01, 0.03, 0.1, 0.3, 1.0) if mode != "none" else (0.0,)):
                    mk = lambda sd, x, v, ns, k=kappa: run_vmn(sd, mode, k, w0, beta, x, v, ns)
                    r = evaluate(mk, TUNE, tune_only=True)
                    grid.append({"kappa": kappa, "val_end_error": float(np.mean([q["val_end_error"] for q in r]))})
                best = min(grid, key=lambda q: q["val_end_error"])
                mk = lambda sd, x, v, ns, k=best["kappa"]: run_vmn(sd, mode, k, w0, beta, x, v, ns)
                r = evaluate(mk, TEST)
                res["vmn"][key] = {"grid": grid, "kappa": best["kappa"],
                                   "test_end_error": float(np.mean([q["test_end_error"] for q in r])),
                                   "test_end_error_sd": float(np.std([q["test_end_error"] for q in r])),
                                   "test_mean_error": float(np.mean([q["test_mean_error"] for q in r])),
                                   "curve": np.mean([q["test_error_curve"] for q in r], 0).tolist()}
                print(key, best["kappa"], round(res["vmn"][key]["test_end_error"], 4),
                      f"{time.time() - t0:.0f}s", flush=True)
    grid = []
    for rho in (0.9, 0.99, 1.05):
        for s in (0.1, 0.5):
            for leak in (0.2, 1.0):
                mk = lambda sd, x, v, ns, a=rho, b=s, c=leak: run_esn(sd, a, b, c, x, v, ns)
                r = evaluate(mk, TUNE, tune_only=True)
                grid.append({"rho": rho, "s": s, "leak": leak,
                             "val_end_error": float(np.mean([q["val_end_error"] for q in r]))})
    best = min(grid, key=lambda q: q["val_end_error"])
    mk = lambda sd, x, v, ns: run_esn(sd, best["rho"], best["s"], best["leak"], x, v, ns)
    r = evaluate(mk, TEST)
    res["esn"] = {"grid": grid, "chosen": best,
                  "test_end_error": float(np.mean([q["test_end_error"] for q in r])),
                  "test_end_error_sd": float(np.std([q["test_end_error"] for q in r])),
                  "test_mean_error": float(np.mean([q["test_mean_error"] for q in r])),
                  "curve": np.mean([q["test_error_curve"] for q in r], 0).tolist()}
    print("esn", best, round(res["esn"]["test_end_error"], 4), flush=True)

    V = res["vmn"]
    verdicts = {}
    for beta in (0.0, 2.0):
        e = lambda m, w: V[f"w0{w:g}_{m}_beta{beta:g}"]["test_end_error"]
        verdicts[f"K1_beta{beta:g}_vortex_beats_linear_and_none"] = bool(e("vortex", 0) < e("linear", 0) and e("vortex", 0) < e("none", 0))
        adv0 = e("none", 0) - e("vortex", 0)
        adv3 = e("none", 3) - e("vortex", 3)
        verdicts[f"K3_beta{beta:g}_advantage_w0=0_vs_w0=3"] = [adv0, adv3]
        verdicts[f"K3_beta{beta:g}_shrinks_by_half"] = bool(adv0 > 0 and adv3 <= 0.5 * adv0)
    best_vmn = min(q["test_end_error"] for q in V.values())
    verdicts["K2_best_vmn_beats_esn"] = bool(best_vmn < res["esn"]["test_end_error"])
    verdicts["best_vmn_end_error"] = best_vmn
    res["verdicts"] = verdicts
    res["seconds"] = round(time.time() - t0, 1)
    os.makedirs("results", exist_ok=True)
    json.dump(res, open("results/gate2_receipt.json", "w"), indent=2)
    print(json.dumps(verdicts, indent=2), res["seconds"])


if __name__ == "__main__":
    main()
