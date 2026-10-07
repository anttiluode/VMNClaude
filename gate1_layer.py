"""Gate 1 - does the vortex structure help a reservoir of oscillator neurons?

Units are Stuart-Landau oscillators (the neuron of VMN_NOTE section 4) at fixed
2D positions p_k, coupled through the vortex strain kernel of section 1:

    dz_k/dt = (mu + i w_k) z_k - (1 + i beta)|z_k|^2 z_k + kappa * C(z)_k + s * win_k * u(t)

    C(z) = M conj(z)   "vortex"  : the antilinear vortex response map (section 1)
    C(z) = M z         "linear"  : same matrix, ordinary linear coupling
    C(z) = 0           "none"    : uncoupled oscillators

    M_kj = i g_j / (2 pi (conj p_k - conj p_j)^2), normalised to spectral norm 1.

Onset is MEASURED, not assumed (Sol): mu_c = -max Re eig of the linearisation at
z = 0 without mu. Every sweep is in delta = mu - mu_c.

Tasks (fixed linear ridge readout on [Re z, Im z, 1]):
  MC     linear memory capacity, sum over delays 1..60 of test R^2
  NARMA  NARMA10 test NRMSE
  CONS   consistency: two copies, same input, different initial state;
         normalised state distance after washout (0 = forgets its start)

Baseline: echo state network with 2N tanh units (same number of real state
variables as N complex oscillators).

Kill conditions, fixed before running (VMN_NOTE section 6, sharpened by Sol):
  K1  vortex coupling must beat linear coupling with the same matrix on MC or
      NARMA; otherwise the vortex structure adds nothing.
  K2  best VMN variant must beat the ESN on MC or NARMA; otherwise not useful.
  K3  the best delta must lie on the stable side close to onset, with
      consistency breaking above onset; otherwise section 4's explanation fails.
  Shear: beta = 0 vs beta = 2 is reported for every variant.
"""
import argparse
import json
import os
import time

import numpy as np

N = 40
TRAIN, VAL, TEST, WASH = 3000, 800, 1200, 200
SUB, DT = 4, 0.25          # 4 RK4 substeps per input, input held for 1 time unit
DELAYS = 60


# ---------------------------------------------------------------------------
# inputs and targets
# ---------------------------------------------------------------------------

def make_inputs(rng, L):
    u_mc = rng.uniform(-0.5, 0.5, L)
    u_na = rng.uniform(0.0, 0.5, L)
    y = np.zeros(L)
    for t in range(9, L - 1):
        y[t + 1] = (0.3 * y[t] + 0.05 * y[t] * y[t - 9:t + 1].sum()
                    + 1.5 * u_na[t - 9] * u_na[t] + 0.1)
    return u_mc, u_na, y


# ---------------------------------------------------------------------------
# VMN reservoir
# ---------------------------------------------------------------------------

def vortex_matrix(rng):
    R = np.sqrt(N / np.pi)
    pts = []
    while len(pts) < N:
        c = rng.uniform(-R, R, 2)
        if c @ c > R * R:
            continue
        if pts and np.min(np.linalg.norm(np.array(pts) - c, axis=1)) < 0.5:
            continue
        pts.append(c)
    p = np.array(pts)
    z = p[:, 0] + 1j * p[:, 1]
    g = rng.choice([-1.0, 1.0], N) * rng.uniform(0.5, 1.0, N)
    d = np.conj(z[:, None] - z[None, :])
    np.fill_diagonal(d, 1.0)
    M = 1j * g[None, :] / (2 * np.pi * d ** 2)
    np.fill_diagonal(M, 0.0)
    return M / np.linalg.norm(M, 2), g


def onset_mu(M, w, kappa, mode):
    """mu_c such that the linearisation at z = 0 is marginal at mu = mu_c."""
    A = np.diag(1j * w)
    if mode == "linear":
        A = A + kappa * M
    # real form of z -> A z + B conj(z)
    B = kappa * M if mode == "vortex" else np.zeros_like(M)
    Rf = np.block([[(A + B).real, (-A + B).imag],
                   [(A + B).imag, (A - B).real]])
    return -float(np.max(np.linalg.eigvals(Rf).real))


def run_vmn(z0, U, M, w, win, mu, beta, kappa, s, mode):
    """z0: (N,B) complex, U: (L,B) inputs. Returns states (L,N,B)."""
    L = U.shape[0]
    lin = (mu + 1j * w)[:, None]
    nl = 1 + 1j * beta
    z = z0.astype(complex)
    out = np.empty((L, N, z.shape[1]), complex)

    def f(z, drive):
        if mode == "vortex":
            c = M @ np.conj(z)
        elif mode == "linear":
            c = M @ z
        else:
            c = 0.0
        return lin * z - nl * (z.real ** 2 + z.imag ** 2) * z + kappa * c + drive

    for t in range(L):
        drive = s * win[:, None] * U[t][None, :]
        for _ in range(SUB):
            k1 = f(z, drive)
            k2 = f(z + 0.5 * DT * k1, drive)
            k3 = f(z + 0.5 * DT * k2, drive)
            k4 = f(z + DT * k3, drive)
            z = z + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if not np.all(np.isfinite(z)):
            out[t:] = np.nan
            return out
        out[t] = z
    return out


def run_esn(x0, U, W, win, s):
    L = U.shape[0]
    x = x0.copy()
    out = np.empty((L, W.shape[0], x.shape[1]))
    for t in range(L):
        x = np.tanh(W @ x + s * win[:, None] * U[t][None, :])
        out[t] = x
    return out


# ---------------------------------------------------------------------------
# readout and scoring
# ---------------------------------------------------------------------------

def features(S):
    """S: (L, n) complex or real -> (L, F)."""
    if np.iscomplexobj(S):
        X = np.concatenate([S.real, S.imag], 1)
    else:
        X = S
    return np.concatenate([X, np.ones((X.shape[0], 1))], 1)


def ridge_fit(X, Y, lam):
    A = X.T @ X + lam * np.eye(X.shape[1])
    return np.linalg.solve(A, X.T @ Y)


LAMS = (1e-8, 1e-6, 1e-4, 1e-2, 1.0)


def score(S_mc, S_na, u_mc, y_na):
    """Train on [WASH, WASH+TRAIN), choose lambda on VAL, report TEST."""
    if not (np.all(np.isfinite(S_mc)) and np.all(np.isfinite(S_na))):
        return {"MC": 0.0, "NARMA": 9.99, "MC_val": 0.0, "NARMA_val": 9.99}
    a, b, c = WASH, WASH + TRAIN, WASH + TRAIN + VAL
    # memory capacity
    X = features(S_mc)
    T = np.stack([np.roll(u_mc, k) for k in range(1, DELAYS + 1)], 1)
    best = None
    for lam in LAMS:
        Wt = ridge_fit(X[a:b], T[a:b], lam)
        pv = X[b:c] @ Wt
        r2v = np.array([max(0.0, np.corrcoef(pv[:, k], T[b:c, k])[0, 1]) ** 2 if np.std(pv[:, k]) > 0 else 0.0
                        for k in range(DELAYS)])
        if best is None or r2v.sum() > best[0]:
            pt = X[c:c + TEST] @ Wt
            r2t = np.array([max(0.0, np.corrcoef(pt[:, k], T[c:c + TEST, k])[0, 1]) ** 2 if np.std(pt[:, k]) > 0 else 0.0
                            for k in range(DELAYS)])
            best = (r2v.sum(), r2t.sum())
    mc_val, mc = best
    # NARMA10
    X = features(S_na)
    best = None
    for lam in LAMS:
        Wt = ridge_fit(X[a:b], y_na[a:b], lam)
        ev = np.sqrt(np.mean((X[b:c] @ Wt - y_na[b:c]) ** 2)) / np.std(y_na[b:c])
        if best is None or ev < best[0]:
            et = np.sqrt(np.mean((X[c:c + TEST] @ Wt - y_na[c:c + TEST]) ** 2)) / np.std(y_na[c:c + TEST])
            best = (ev, et)
    na_val, na = best
    return {"MC": float(mc), "NARMA": float(na), "MC_val": float(mc_val), "NARMA_val": float(na_val)}


def consistency(Sa, Sb):
    if not (np.all(np.isfinite(Sa)) and np.all(np.isfinite(Sb))):
        return float("nan")
    tail = slice(WASH + TRAIN, None)
    return float(np.linalg.norm(Sa[tail] - Sb[tail]) / max(np.linalg.norm(Sa[tail]), 1e-12))


# ---------------------------------------------------------------------------
# one system, one configuration -> all three measurements
# ---------------------------------------------------------------------------

def eval_vmn(seed, mode, beta, kappa, s, delta, counter=False):
    """counter=True: each unit rotates in the sense of its circulation, w_k -> sign(g_k) w_k
    (mixed co- and counter-rotating units, as with vortices of both signs)."""
    rng = np.random.default_rng(1000 + seed)
    M, g = vortex_matrix(rng)
    w = rng.uniform(0.5, 1.5, N)
    if counter:
        w = w * np.sign(g)
    win = rng.normal(size=N) + 1j * rng.normal(size=N)
    win /= np.abs(win)
    L = WASH + TRAIN + VAL + TEST
    u_mc, u_na, y_na = make_inputs(np.random.default_rng(seed), L)
    mu_c = onset_mu(M, w, kappa if mode != "none" else 0.0, mode)
    mu = mu_c + delta
    U = np.stack([u_mc, u_na, u_mc], 1)
    z0 = np.zeros((N, 3), complex)
    z0[:, 2] = 0.3 * (rng.normal(size=N) + 1j * rng.normal(size=N))   # different start
    z0[:, 0] = 0.3 * (rng.normal(size=N) + 1j * rng.normal(size=N))
    S = run_vmn(z0, U, M, w, win, mu, beta, kappa, s, mode)
    r = score(S[:, :, 0], S[:, :, 1], u_mc, y_na)
    r["CONS"] = consistency(S[:, :, 0], S[:, :, 2])
    r["mu_c"] = mu_c
    return r


def eval_esn(seed, rho, s):
    rng = np.random.default_rng(2000 + seed)
    n = 2 * N
    W = rng.normal(size=(n, n))
    W *= rho / np.max(np.abs(np.linalg.eigvals(W)))
    win = rng.uniform(-1, 1, n)
    L = WASH + TRAIN + VAL + TEST
    u_mc, u_na, y_na = make_inputs(np.random.default_rng(seed), L)
    U = np.stack([u_mc, u_na, u_mc], 1)
    x0 = np.zeros((n, 3))
    x0[:, 2] = rng.uniform(-1, 1, n)
    x0[:, 0] = rng.uniform(-1, 1, n)
    S = run_esn(x0, U, W, win, s)
    r = score(S[:, :, 0], S[:, :, 1], u_mc, y_na)
    r["CONS"] = consistency(S[:, :, 0], S[:, :, 2])
    return r


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    tune_seeds = (0,) if args.quick else (0, 1)
    sweep_seeds = (10, 11) if args.quick else (10, 11, 12)
    deltas = (-0.3, -0.1, -0.03, -0.01, 0.0, 0.03, 0.1) if args.quick else \
             (-0.5, -0.2, -0.1, -0.05, -0.02, -0.01, 0.0, 0.01, 0.03, 0.1)
    variants = [(m, b) for m in ("vortex", "linear", "none") for b in (0.0, 2.0)]
    grid_k = (0.3, 1.0, 3.0)
    grid_s = (0.1, 0.3, 1.0)
    out = {"setup": {"N_oscillators": N, "real_state_dim": 2 * N, "train": TRAIN, "val": VAL,
                     "test": TEST, "washout": WASH, "substeps": SUB, "dt": DT, "delays": DELAYS,
                     "tune_seeds": tune_seeds, "sweep_seeds": sweep_seeds,
                     "tuning_delta": -0.05, "deltas": deltas},
           "tuning": {}, "sweep": {}, "esn": {}}

    # 1) tune kappa and input scale per variant and per task, on validation only
    chosen = {}
    for mode, beta in variants:
        key = f"{mode}_beta{beta:g}"
        rows = []
        for k in (grid_k if mode != "none" else (0.0,)):
            for s in grid_s:
                res = [eval_vmn(sd, mode, beta, k, s, -0.05) for sd in tune_seeds]
                rows.append({"kappa": k, "s": s,
                             "MC_val": float(np.mean([r["MC_val"] for r in res])),
                             "NARMA_val": float(np.mean([r["NARMA_val"] for r in res]))})
        bmc = max(rows, key=lambda r: r["MC_val"])
        bna = min(rows, key=lambda r: r["NARMA_val"])
        chosen[key] = {"MC": (bmc["kappa"], bmc["s"]), "NARMA": (bna["kappa"], bna["s"])}
        out["tuning"][key] = {"grid": rows, "chosen": chosen[key]}
        print(key, "tuned", chosen[key], f"{time.time() - t0:.0f}s", flush=True)

    # 2) sweep distance to the measured onset with the tuned settings
    for mode, beta in variants:
        key = f"{mode}_beta{beta:g}"
        rows = []
        for d in deltas:
            row = {"delta": d}
            for task in ("MC", "NARMA"):
                k, s = chosen[key][task]
                res = [eval_vmn(sd, mode, beta, k, s, d) for sd in sweep_seeds]
                row[task] = float(np.mean([r[task] for r in res]))
                row[task + "_sd"] = float(np.std([r[task] for r in res]))
                row["CONS_" + task + "_cfg"] = float(np.nanmean([r["CONS"] for r in res]))
                row["mu_c_" + task + "_cfg"] = float(np.mean([r["mu_c"] for r in res]))
            rows.append(row)
        out["sweep"][key] = rows
        print(key, "swept", f"{time.time() - t0:.0f}s", flush=True)

    # 3) ESN baseline, tuned the same way
    rows = []
    for rho in (0.5, 0.9, 1.1):
        for s in (0.1, 0.5, 1.0):
            res = [eval_esn(sd, rho, s) for sd in tune_seeds]
            rows.append({"rho": rho, "s": s,
                         "MC_val": float(np.mean([r["MC_val"] for r in res])),
                         "NARMA_val": float(np.mean([r["NARMA_val"] for r in res]))})
    bmc = max(rows, key=lambda r: r["MC_val"])
    bna = min(rows, key=lambda r: r["NARMA_val"])
    esn = {}
    for task, b in (("MC", bmc), ("NARMA", bna)):
        res = [eval_esn(sd, b["rho"], b["s"]) for sd in sweep_seeds]
        esn[task] = float(np.mean([r[task] for r in res]))
        esn[task + "_sd"] = float(np.std([r[task] for r in res]))
        esn[task + "_cfg"] = {"rho": b["rho"], "s": b["s"]}
        esn["CONS_" + task + "_cfg"] = float(np.mean([r["CONS"] for r in res]))
    out["esn"] = {"grid": rows, "test": esn}

    # 4) verdicts
    def best_of(key, task):
        rows = out["sweep"][key]
        return (max if task == "MC" else min)(rows, key=lambda r: r[task])

    summ = {}
    for mode, beta in variants:
        key = f"{mode}_beta{beta:g}"
        summ[key] = {"MC_best": best_of(key, "MC")["MC"], "MC_best_delta": best_of(key, "MC")["delta"],
                     "NARMA_best": best_of(key, "NARMA")["NARMA"], "NARMA_best_delta": best_of(key, "NARMA")["delta"]}
    out["summary"] = summ
    out["summary"]["esn"] = {"MC": esn["MC"], "NARMA": esn["NARMA"]}
    v = {}
    for beta in (0.0, 2.0):
        a, b = summ[f"vortex_beta{beta:g}"], summ[f"linear_beta{beta:g}"]
        v[f"K1_beta{beta:g}_vortex_beats_linear"] = {
            "MC": a["MC_best"] > b["MC_best"], "NARMA": a["NARMA_best"] < b["NARMA_best"]}
    vm = [summ[k] for k in summ if k != "esn"]
    v["K2_best_vmn_beats_esn"] = {"MC": max(r["MC_best"] for r in vm) > esn["MC"],
                                  "NARMA": min(r["NARMA_best"] for r in vm) < esn["NARMA"]}
    v["K3_best_delta_by_variant"] = {k: (summ[k]["MC_best_delta"], summ[k]["NARMA_best_delta"])
                                    for k in summ if k != "esn"}
    out["verdicts"] = v
    out["seconds"] = round(time.time() - t0, 1)
    os.makedirs("results", exist_ok=True)
    with open("results/gate1_receipt.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({"summary": out["summary"], "verdicts": v, "seconds": out["seconds"]}, indent=2))


if __name__ == "__main__":
    main()
