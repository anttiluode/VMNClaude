"""Gate 6 - read the goal, listen, then restore the memory.

Set up after Sol's review of Gate 5 (VMN repo). Sol's three conditions for a useful result:
  1. measure retained phase disturbance against a same-time UNQUERIED control, alongside
     position error;
  2. compare compensating schemes, counting every pulse and the full listening time;
  3. show the compensation keeps the useful query answer, using observations from ONE bank.
And: "merely making disturbance smaller can be achieved by making the query weaker", so
schemes are compared at MATCHED ANSWER ACCURACY, across a sweep of ping amplitudes.

Store: path integration exactly as Gates 2 and 4 (40 uncoupled Stuart-Landau units, beta = 0,
sigma = 0.05, 200 steps). Query: a goal q = x + D, |D| <= 0.35, arrives as a kick in the phase
pattern the bank would hold at q (Gate 4).

Schemes (each window = 16 substeps = 4 time units of listening, standing still):
  single        +a e^{i phi}                                  1 pulse,  1 window
  repeat        +a e^{i phi}, then +a e^{i phi} again          2 pulses, 2 windows (no compensation;
                                                                controls for "two looks average noise")
  sign_flip     +a e^{i phi}, then -a e^{i phi}                2 pulses, 2 windows (Sol's compensator;
                                                                residual ~ a^2 sin 2 alpha)
  conjugate     +a e^{i phi}, then a u^2 e^{-i phi}            2 pulses, 2 windows. u = z/|z|: the second
                                                                kick sits at relative phase -alpha, so with
                                                                full relaxation it cancels the first EXACTLY
                                                                (arg(1+a e^{ia}) is odd). Physically this is
                                                                a phase-conjugating (2f, parametric) drive;
                                                                the kick depends on the bank's own state.
Query code: "uniform" (one D for all units, as Gate 4). A "scaled" code (each module pinged with
D * period_m / period_coarse, as grid sweeps scale with module spacing, Vollan et al. 2025) was
written and dropped BEFORE running: building that pattern for a goal given in world coordinates
needs the true position, which hands the answer to the reader and bypasses the memory. Grid sweeps
are self-referenced (Gate 5's setting), not world-coordinate goal queries.

ONE-BANK observations (no noise-cancelling twin): before the ping the bank's phases theta_0 are
read once; after each pulse, at t = 0.25, 0.5, 1, 2, 4: amplitude deviation |z| - sqrt(mu) and
wrapped phase change theta(t) - theta_0. Readers: ridge and kNN (Gate 4), best on validation.

Damage (measurement only): a twin that is never pinged shares the noise draws in each copy's own
phase frame (Gate 5's co-rotating noise). After the last window plus 16 substeps of relaxation:
RMS wrapped per-unit phase difference, and the least-squares position shift.

Amplitudes a in {0.03, 0.05, 0.1, 0.2, 0.3, 0.5}. For each scheme this gives a curve of
(answer error, damage). Damage of each pair scheme is compared with `single` at the same answer
error (log-log interpolation along single's curve).

Kill conditions, fixed before running:
  K1  restoration is real: at the answer error `single` reaches with a = 0.3, sign_flip OR
      conjugate retains <= 0.5x the phase disturbance single needs for that error, in the query code. Otherwise read-and-restore buys nothing over a weaker query.
  K2  it is the compensation, not the second look: the winning pair scheme also retains <= 0.5x
      the disturbance of `repeat` at matched answer error.
  K3  exactness matters: conjugate beats sign_flip by >= 2x at matched answer error.
"""
import json
import os
import time

import numpy as np

import gate2_path as g2
import gate4_ping as g4
import gate5_sweeps as g5

N, MU, SIGMA, DT, SUB = g2.N, g2.MU, g2.SIGMA, g2.DT, g2.SUB
E = {"tr": 1000, "va": 300, "te": 500}
OFF = {"tr": 100, "va": 200, "te": 300}
SEEDS = (10, 11, 12)
AMPS = (0.03, 0.05, 0.1, 0.2, 0.3, 0.5)
SCHEMES = ("single", "repeat", "sign_flip", "conjugate")
CODES = ("uniform",)   # "scaled" dropped before running: see module docstring
H = 16
SAMPLE = (1, 2, 4, 8, 16)
RELAX = 16
PERIOD = {2 * np.pi * 1.0: 1.0, 2 * np.pi * 2.0: 0.5, 2 * np.pi * 4.0: 0.25}


def step(zp, zu, rng):
    ep, eu = g5.shared(g5.noise(rng, zp.shape), zp, zu)
    return g5.rk4(zp, MU, 0.0) + ep, g5.rk4(zu, MU, 0.0) + eu


def query_phase(code, d, b, q, x_hat):
    """phi_k: phase pattern of the goal. scaled: each module sees D shrunk to its own period."""
    if code == "uniform":
        return b[:, None] * (d @ q.T)
    scale = np.array([PERIOD[bk] for bk in b]) / 1.0                  # coarse period = 1
    D = q - x_hat                                                       # true offset (task input)
    return b[:, None] * (d @ x_hat.T) + b[:, None] * scale[:, None] * (d @ D.T)


def run(sd, Z, q, x, code, scheme, amp, rng):
    _, d, b = g2.bank(sd)
    # the query pattern is defined in world coordinates for "uniform"; for "scaled" the module
    # offsets are relative to the true position x (the brain sweeps from its own map position;
    # here the map error is small next to the scaled offsets, and x only sets the pattern)
    phi = query_phase(code, d, b, q, x)
    zp, zu = Z.copy(), Z.copy()
    th0 = np.angle(zp)
    feats = []
    pulses = [amp * np.sqrt(MU) * np.exp(1j * phi)]
    if scheme != "single":
        pulses.append(None)
    for k, _ in enumerate(pulses):
        if k == 0:
            kick = amp * np.sqrt(MU) * np.exp(1j * phi)
        else:
            u = zp / np.maximum(np.abs(zp), 1e-12)
            if scheme == "repeat":
                kick = amp * np.sqrt(MU) * np.exp(1j * phi)
            elif scheme == "sign_flip":
                kick = -amp * np.sqrt(MU) * np.exp(1j * phi)
            else:                                                  # conjugate
                kick = amp * np.sqrt(MU) * u * u * np.exp(-1j * phi)
        zp = zp + kick
        for s in range(1, H + 1):
            zp, zu = step(zp, zu, rng)
            if s in SAMPLE:
                feats.append(np.abs(zp) - np.sqrt(MU))
                feats.append(np.angle(zp * np.exp(-1j * th0)))
    for _ in range(RELAX):
        zp, zu = step(zp, zu, rng)
    dphi = np.angle(zp * np.conj(zu))
    shift = np.linalg.lstsq(b[:, None] * d, dphi, rcond=None)[0].T
    return (np.concatenate(feats, 0).T,
            float(np.sqrt((dphi ** 2).mean())),
            float(np.linalg.norm(shift, axis=1).mean()))


def run_seed(sd):
    data, Z = {}, {}
    for s in E:
        xx, v = g2.trajectories(np.random.default_rng(OFF[s] + sd), E[s])
        rng = np.random.default_rng(900 + OFF[s] + sd)
        ang = rng.uniform(0, 2 * np.pi, E[s]); rad = 0.35 * np.sqrt(rng.uniform(0, 1, E[s]))
        D = np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
        data[s] = {"x": xx[-1], "D": D, "q": xx[-1] + D}
        Z[s] = g4.integrate(sd, "none", 0.0, xx, v, np.random.default_rng(31 + 7 * sd + {"tr": 1, "va": 2, "te": 3}[s]))
    Y = {s: data[s]["D"] for s in E}
    res = {"predict_zero_error": float(np.linalg.norm(Y["te"], axis=1).mean())}
    Fo = {s: g4.oracle_feats(sd, Z[s], data[s]["q"]) for s in E}
    res["oracle_uniform"] = min(g4.score(Fo, Y).values())
    for code in CODES:
        for scheme in SCHEMES:
            for amp in AMPS:
                F, dmg, sh = {}, None, None
                for s in E:
                    rng = np.random.default_rng(5000 + 13 * sd + {"tr": 1, "va": 2, "te": 3}[s])
                    F[s], ph, shf = run(sd, Z[s], data[s]["q"], data[s]["x"], code, scheme, amp, rng)
                    if s == "te":
                        dmg, sh = ph, shf
                sc = g4.score(F, Y)
                res[f"{code}|{scheme}|{amp:g}"] = {"answer_error": min(sc.values()), **{f"err_{k}": v for k, v in sc.items()},
                                                   "phase_disturbance_rms": dmg, "position_shift": sh}
    return res


def interp_damage(curve, target):
    """Damage single would need to reach answer error `target`, interpolating log-log along its curve.
    curve: list of (answer_error, damage), any order. Returns None if target outside its range."""
    pts = sorted(curve)
    e = np.log([p[0] for p in pts]); dmg = np.log([p[1] for p in pts])
    if not (e.min() <= np.log(target) <= e.max()):
        return None
    order = np.argsort(e)
    return float(np.exp(np.interp(np.log(target), e[order], dmg[order])))


def main():
    t0 = time.time()
    per = {}
    for sd in SEEDS:
        per[sd] = run_seed(sd)
        print("seed", sd, f"{time.time() - t0:.0f}s", flush=True)
    keys = [k for k in per[SEEDS[0]] if "|" in k]
    S = {"predict_zero_error": float(np.mean([per[s]["predict_zero_error"] for s in SEEDS])),
         "oracle_uniform": float(np.mean([per[s]["oracle_uniform"] for s in SEEDS]))}
    for k in keys:
        S[k] = {m: [float(np.mean([per[s][k][m] for s in SEEDS])), float(np.std([per[s][k][m] for s in SEEDS]))]
                for m in per[SEEDS[0]][k]}

    matched = {}
    for code in CODES:
        curves = {sch: [(S[f"{code}|{sch}|{a:g}"]["answer_error"][0], S[f"{code}|{sch}|{a:g}"]["phase_disturbance_rms"][0])
                        for a in AMPS] for sch in SCHEMES}
        target = S[f"{code}|single|0.3"]["answer_error"][0]
        row = {"target_answer_error": target}
        for sch in SCHEMES:
            row[sch] = interp_damage(curves[sch], target)
        matched[code] = row
        print(code, "curves:")
        for sch in SCHEMES:
            print("  ", sch.ljust(10), " ".join(f"({e:.4f},{dmg:.4f})" for e, dmg in curves[sch]))
        print("  matched at", round(target, 4), {k: (round(v, 5) if v else v) for k, v in row.items()})

    def r(code, a, b):
        x, y = matched[code][a], matched[code][b]
        return None if (x is None or y is None) else y / x

    best = {}
    for code in CODES:
        cands = {s: matched[code][s] for s in ("sign_flip", "conjugate") if matched[code][s] is not None}
        best[code] = min(cands, key=cands.get) if cands else None
    V = {"K1_restore_beats_weaker_query": bool(all(best[c] and r(c, best[c], "single") >= 2 for c in CODES)),
         "K1_single_over_best": {c: [best[c], r(c, best[c], "single") if best[c] else None] for c in CODES},
         "K2_compensation_not_second_look": bool(all(best[c] and r(c, best[c], "repeat") and r(c, best[c], "repeat") >= 2 for c in CODES)),
         "K2_repeat_over_best": {c: r(c, best[c], "repeat") if best[c] else None for c in CODES},
         "K3_conjugate_beats_sign_flip_2x": bool(all(r(c, "conjugate", "sign_flip") and r(c, "conjugate", "sign_flip") >= 2 for c in CODES)),
         "K3_signflip_over_conjugate": {c: r(c, "conjugate", "sign_flip") for c in CODES}}
    out = {"setup": {"N": N, "episodes": E, "seeds": SEEDS, "amps": AMPS, "schemes": SCHEMES, "codes": CODES,
                     "listen_substeps_per_pulse": H, "sample_substeps": SAMPLE, "relax_substeps": RELAX},
           "summary_mean_sd": S, "matched": matched, "verdicts": V,
           "per_seed": {str(k): v for k, v in per.items()}, "seconds": round(time.time() - t0, 1)}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/gate6_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), out["seconds"])


if __name__ == "__main__":
    main()
