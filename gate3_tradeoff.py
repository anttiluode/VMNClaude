"""Gate 3 - the memory / control trade-off, and event-gated coupling.

Claim under test (from the Oct 7 discussion):
  A phase coordinate cannot be both a perfect memory and a knob on the response operator.
  An uncoupled bank has symmetry U(1)^N (rotate any unit's phase alone): every phase is a
  neutral direction, so it integrates its input perfectly, but the response operator at two
  stored positions differs only by per-unit rotations, so its singular values cannot depend on
  the stored position. Coupling breaks the symmetry: the operator starts to "see" the stored
  position (control) and the phases stop being clean integrators (drift).
  Refinement: ordinary coupling Mz still keeps the GLOBAL U(1) (rotate all phases together);
  vortex coupling M conj(z) breaks it too. Path integration stores position in per-unit phase
  differences, so here both couplings are expected to trade.

Two measured axes, swept over kappa for vortex and linear coupling (beta = 0, w0 = 0):
  ERROR   - Gate 2 end-of-episode position error (noise sigma = 0.05, readout ridge on phases).
  CONTROL - how much the stored position changes the response operator. For K = 48 test
            positions, put every unit on its limit cycle at the true phase b_k d_k.x, form the
            exact tangent D (doubled complex form, A/B blocks), take Phi = expm(D H), H = 2,
            and its singular values s(x). CONTROL = mean pairwise ||s(x)-s(x')|| / mean ||s||.
            Rotation-invariant: at kappa = 0 it is exactly 0 by symmetry.

Event-gated variant: kappa = 0 between events, kappa = kappa_on for W = 5 steps every 50
steps (duty 10%). Its CONTROL during an event is that of kappa_on; its ERROR is measured.

Kill conditions, fixed before running:
  K1  trade-off: CONTROL(kappa=0) < 1e-9, CONTROL rises with kappa, and no kappa > 0 has
      ERROR below ERROR(kappa=0) by more than 2 seed-sd. A coupling that both helps
      integration and gives control refutes the claim.
  K2  global-phase refinement: rotating all phases together leaves s unchanged (< 1e-9 rel)
      for linear coupling and changes it (> 1e-3 rel) for vortex coupling.
  K3  event gating: for at least one kappa_on, in BOTH modes, the gated variant's CONTROL is
      >= 2x the CONTROL of fixed coupling at the same ERROR (fixed curve interpolated in
      log kappa). Otherwise gating buys nothing over simply using a weaker coupling.
"""
import json
import os
import time

import numpy as np
from scipy.linalg import expm

import gate2_path as g2

N, MU, SIGMA = g2.N, g2.MU, g2.SIGMA
KAPPAS = (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)
KON = (0.03, 0.1, 0.3, 1.0)
EVERY, WIN = 50, 5
H = 2.0
TEST = (10, 11, 12)


# ---------- control axis: rotation-invariant operator change -------------------------------
def tangent(z, M, kappa, mode, beta=0.0):
    """Exact tangent in doubled form: d(dz)/dt = A dz + B conj(dz)."""
    c = 1 + 1j * beta
    A = np.diag(MU - 2 * c * np.abs(z) ** 2).astype(complex)
    B = np.diag(-c * z ** 2)
    if mode == "linear":
        A = A + kappa * M
    elif mode == "vortex":
        B = B + kappa * M
    return np.block([[A, B], [np.conj(B), np.conj(A)]])


def svals(z, M, kappa, mode):
    return np.linalg.svd(expm(tangent(z, M, kappa, mode) * H), compute_uv=False)


def control(seed, kappa, mode, K=48):
    M, d, b = g2.bank(seed)
    xs = np.random.default_rng(500 + seed).uniform(0.1, 0.9, (K, 2))
    S = np.array([svals(np.sqrt(MU) * np.exp(1j * b * (d @ x)), M, kappa, mode) for x in xs])
    diffs = np.linalg.norm(S[:, None] - S[None], axis=-1)[np.triu_indices(K, 1)]
    return float(diffs.mean() / np.linalg.norm(S, axis=1).mean())


def global_rotation_check(seed=10, kappa=0.3):
    M, d, b = g2.bank(seed)
    x = np.array([0.4, 0.6])
    z = np.sqrt(MU) * np.exp(1j * b * (d @ x))
    out = {}
    for mode in ("linear", "vortex"):
        s0 = svals(z, M, kappa, mode)
        rel = max(np.linalg.norm(svals(z * np.exp(1j * th), M, kappa, mode) - s0) / np.linalg.norm(s0)
                  for th in np.linspace(0.3, np.pi, 6))
        out[mode] = float(rel)
    return out


# ---------- error axis: path integration, fixed or event-gated coupling --------------------
def run_gated(seed, mode, kappa_on, x, v, noise_seed):
    """Gate 2 dynamics (beta = 0, w0 = 0) with kappa on only inside event windows."""
    M, d, b = g2.bank(seed)
    E = x.shape[1]
    rng = np.random.default_rng(noise_seed)
    z = np.sqrt(MU) * np.exp(1j * b[:, None] * (d @ x[0].T))
    feats = np.empty((g2.STEPS + 1, E, 2 * N))

    def feat(z):
        ph = z / np.maximum(np.abs(z), 1e-9)
        return np.concatenate([ph.real, ph.imag], 0).T

    feats[0] = feat(z)
    DT, SUB = g2.DT, g2.SUB
    for t in range(g2.STEPS):
        kappa = kappa_on if (t % EVERY) >= EVERY - WIN else 0.0
        lin = MU + 1j * (b[:, None] * (d @ v[t].T) / (SUB * DT))

        def f(z):
            c = M @ np.conj(z) if mode == "vortex" else M @ z
            return lin * z - (z.real ** 2 + z.imag ** 2) * z + kappa * c

        for _ in range(SUB):
            k1 = f(z); k2 = f(z + 0.5 * DT * k1)
            k3 = f(z + 0.5 * DT * k2); k4 = f(z + DT * k3)
            z = z + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            z = z + SIGMA * np.sqrt(DT) * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape)) / np.sqrt(2)
        feats[t + 1] = feat(z)
    return feats


def error(make):
    r = g2.evaluate(make, TEST)
    e = [q["test_end_error"] for q in r]
    return float(np.mean(e)), float(np.std(e))


def interp_fixed_control(curve, err):
    """CONTROL of fixed coupling at a given ERROR, interpolated along the (error, control) curve."""
    pts = sorted((c["error"], c["control"]) for c in curve)
    es, cs = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    if err <= es[0]:
        return float(cs[0])
    if err >= es[-1]:
        return float(cs[-1])
    return float(np.interp(err, es, cs))


def main():
    t0 = time.time()
    res = {"setup": {"kappas": KAPPAS, "kappa_on": KON, "event_every": EVERY, "window": WIN,
                     "horizon_H": H, "sigma": SIGMA, "test_seeds": TEST,
                     "error_if_always_predicting_box_centre": 0.3826},
           "global_rotation_rel_change": global_rotation_check(), "fixed": {}, "gated": {}}
    print("global rotation:", res["global_rotation_rel_change"], flush=True)

    for mode in ("vortex", "linear"):
        curve = []
        for k in KAPPAS:
            mk = lambda sd, x, v, ns, k=k: g2.run_vmn(sd, mode if k > 0 else "none", k, 0.0, 0.0, x, v, ns)
            e, s = error(mk)
            c = float(np.mean([control(sd, k, mode) for sd in TEST]))
            curve.append({"kappa": k, "error": e, "error_sd": s, "control": c})
            print(mode, "fixed", k, round(e, 4), round(c, 5), f"{time.time() - t0:.0f}s", flush=True)
        res["fixed"][mode] = curve
        gated = []
        for k in KON:
            mk = lambda sd, x, v, ns, k=k: run_gated(sd, mode, k, x, v, ns)
            e, s = error(mk)
            c = next(q["control"] for q in curve if q["kappa"] == k)
            cf = interp_fixed_control(curve, e)
            gated.append({"kappa_on": k, "error": e, "error_sd": s, "control": c,
                          "fixed_control_at_same_error": cf,
                          "ratio": float(c / cf) if cf > 0 else float("inf")})
            print(mode, "gated", k, round(e, 4), round(c, 5), "fixed@err", round(cf, 5), f"{time.time() - t0:.0f}s", flush=True)
        res["gated"][mode] = gated

    V = {}
    for mode in ("vortex", "linear"):
        cur = res["fixed"][mode]
        base = cur[0]
        ctrl = [q["control"] for q in cur]
        V[f"K1_{mode}_control_zero_at_kappa0"] = bool(ctrl[0] < 1e-9)
        V[f"K1_{mode}_control_rises"] = bool(all(b2 >= a2 for a2, b2 in zip(ctrl[1:], ctrl[2:])) and ctrl[1] > ctrl[0])
        V[f"K1_{mode}_no_coupling_beats_kappa0"] = bool(not any(
            q["error"] < base["error"] - 2 * max(base["error_sd"], q["error_sd"]) for q in cur[1:]))
        V[f"K1_{mode}_pass"] = bool(V[f"K1_{mode}_control_zero_at_kappa0"] and V[f"K1_{mode}_control_rises"]
                                    and V[f"K1_{mode}_no_coupling_beats_kappa0"])
        V[f"K3_{mode}_best_ratio"] = max(q["ratio"] for q in res["gated"][mode])
    g = res["global_rotation_rel_change"]
    V["K2_pass"] = bool(g["linear"] < 1e-9 and g["vortex"] > 1e-3)
    V["K3_pass"] = bool(V["K3_vortex_best_ratio"] >= 2 and V["K3_linear_best_ratio"] >= 2)
    res["verdicts"] = V
    res["seconds"] = round(time.time() - t0, 1)
    os.makedirs("results", exist_ok=True)
    json.dump(res, open("results/gate3_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), res["seconds"])


if __name__ == "__main__":
    main()
