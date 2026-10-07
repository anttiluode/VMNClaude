"""Gate 5 - does left-right alternation protect a phase memory from its own reads?

Hypothesis (Claude, 7 Oct 2026, from Sol's PING.md + Vollan et al. Nature 2025):
  Sol measured that pinging an oscillator bank costs state, and the cost accumulates
  (one ping +0.0016 position error, four pings +0.0084). Grid cells "ping" their own map
  every theta cycle with sweeps that alternate ~30 deg left / ~30 deg right of heading.
  Maybe the alternation exists to cancel the read cost, the way spin echo cancels dephasing.

Model: Gate 2's bank (40 Stuart-Landau units, three spatial scales, sigma = 0.05).
A sweep to offset D is a SELF-REFERENCED ping: the kick is a*sqrt(mu)*u_k*exp(i b_k d_k.D),
where u_k = z_k/|z_k| is the unit's current phase. So the sweep goes out from wherever the
map currently thinks it is, as a brain sweep would. A pinged and an unpinged copy share
noise, so their phase difference is exactly the damage the reads did.

Damage is measured as a POSITION SHIFT: least-squares fit of the wrapped per-unit phase
differences to b_k d_k . dx (the shift of the decoded position). Also split into the
component along heading (forward) and across it (lateral).

Exact fact the gate leans on (beta = 0, amplitude fully relaxed between pings):
  a ping at relative phase alpha moves the unit's phase by f(alpha) = arg(1 + a e^{i alpha}),
  an ODD function of alpha. A ping at +D followed by one at -D therefore cancels exactly,
  per unit, to all orders. Left/right sweeps (+30, -30 deg about heading) are NOT negatives
  of each other: their lateral parts are mirror images but their forward parts are equal.
  With shear (beta != 0) the amplitude kick a cos(alpha) also turns into phase; that part is
  EVEN in alpha, so even antipodal pairs stop cancelling.

Sweep patterns (length s = 0.2, heading h):
  same        every ping at R(+30) h s
  left_right  alternating R(+30) h s, R(-30) h s        (the brain's pattern)
  fwd_back    alternating +h s, -h s
  antipodal   alternating +R(30) h s, -R(30) h s
  random      uniform random direction each ping

Part A (stationary): bank held at a random place, n pings (n = 1..32), interval 8 or 2
  substeps between pings, 16 substeps to relax, beta in {0, 0.5}, a in {0.1, 0.3}.
Part B (navigating): Gate 2 trajectories, one sweep every 2 steps along the current
  velocity heading (100 sweeps per episode), beta = 0, a in {0.1, 0.3}. Damage = final
  pinged-vs-unpinged position shift.

Kill conditions, fixed before running:
  H1 (my hypothesis)  left_right cuts mean |shift| vs same by >= 2x, at n = 16 stationary
                      (beta 0, interval 8, a = 0.3) AND in Part B at a = 0.3. Otherwise dead.
  H2 (theory check)   antipodal cuts it vs same by >= 10x in the same two settings.
  H3 (shear)          with beta = 0.5, antipodal's advantage over same falls below 10x.
"""
import json
import os
import time

import numpy as np

import gate2_path as g2

N, MU, SIGMA, DT, SUB, STEPS = g2.N, g2.MU, g2.SIGMA, g2.DT, g2.SUB, g2.STEPS
S_LEN = 0.2
ANGLE = np.deg2rad(30)
SEEDS = (10, 11, 12)
E_A, E_B = 400, 300
NS = (1, 2, 4, 8, 16, 32)
PATTERNS = ("none", "same", "left_right", "fwd_back", "antipodal", "random")
NOISE_FRAME = "corotating"   # "cartesian" reproduces the first run (results/gate5_receipt_cartesian_noise.json)


def rot(h, ang):
    c, s = np.cos(ang), np.sin(ang)
    return np.stack([c * h[:, 0] - s * h[:, 1], s * h[:, 0] + c * h[:, 1]], 1)


def offsets(pattern, k, h, rng):
    """Sweep offset (E,2) for the k-th ping (k = 0, 1, ...) given heading h (E,2)."""
    if pattern == "none":
        return None
    if pattern == "same":
        u = rot(h, ANGLE)
    elif pattern == "left_right":
        u = rot(h, ANGLE if k % 2 == 0 else -ANGLE)
    elif pattern == "fwd_back":
        u = h if k % 2 == 0 else -h
    elif pattern == "antipodal":
        u = rot(h, ANGLE) * (1 if k % 2 == 0 else -1)
    elif pattern == "random":
        a = rng.uniform(0, 2 * np.pi, len(h))
        u = np.stack([np.cos(a), np.sin(a)], 1)
    else:
        raise ValueError(pattern)
    return S_LEN * u


def rhs(z, lin, beta):
    return lin * z - (1 + 1j * beta) * (z.real ** 2 + z.imag ** 2) * z


def rk4(z, lin, beta):
    k1 = rhs(z, lin, beta); k2 = rhs(z + 0.5 * DT * k1, lin, beta)
    k3 = rhs(z + 0.5 * DT * k2, lin, beta); k4 = rhs(z + DT * k3, lin, beta)
    return z + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def noise(rng, shape):
    return SIGMA * np.sqrt(DT) * (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / np.sqrt(2)


def shared(e, zp, zu):
    """Apply one noise draw to both copies. Co-rotating: the draw is expressed in each copy's
    own phase frame (radial, tangential), so a copy at a different phase gets the SAME phase
    kick. Isotropic complex noise rotated by any phase has the same distribution, so this is
    the same noise process; it just stops counting noise reshuffling as read damage."""
    if NOISE_FRAME == "cartesian":
        return e, e
    up = zp / np.maximum(np.abs(zp), 1e-12)
    uu = zu / np.maximum(np.abs(zu), 1e-12)
    return e * up, e * uu


def ping(z, D, amp, d, b):
    if D is None:
        return z
    u = z / np.maximum(np.abs(z), 1e-12)
    return z + amp * np.sqrt(MU) * u * np.exp(1j * b[:, None] * (d @ D.T))


def shift_of(zp, zu, d, b):
    """Least-squares position shift explaining the per-unit phase differences. (E,2)"""
    dphi = np.angle(zp * np.conj(zu))                       # (N,E) wrapped
    A = b[:, None] * d                                      # (N,2)
    return np.linalg.lstsq(A, dphi, rcond=None)[0].T


def split(shift, h):
    fwd = (shift * h).sum(1)
    lat = shift[:, 0] * h[:, 1] * -1 + shift[:, 1] * h[:, 0]
    return fwd, lat


def stats(shift, h):
    fwd, lat = split(shift, h)
    return {"mean_abs_shift": float(np.linalg.norm(shift, axis=1).mean()),
            "mean_forward": float(fwd.mean()), "mean_abs_lateral": float(np.abs(lat).mean()),
            "mean_lateral": float(lat.mean())}


# ---------- Part A: stationary --------------------------------------------------------------
def part_a(sd, beta, amp, interval):
    _, d, b = g2.bank(sd)
    rng0 = np.random.default_rng(500 + sd)
    x = rng0.uniform(0.2, 0.8, (E_A, 2))
    a = rng0.uniform(0, 2 * np.pi, E_A)
    h = np.stack([np.cos(a), np.sin(a)], 1)
    out = {}
    for pat in PATTERNS:
        z0 = np.sqrt(MU) * np.exp(1j * b[:, None] * (d @ x.T))
        # resting state under shear turns at -beta*mu; both copies share it, so diffs are fine
        zp, zu = z0.copy(), z0.copy()
        nrng = np.random.default_rng(40 + sd)
        prng = np.random.default_rng(60 + sd)
        res, k = {}, 0
        for n_target in NS:
            while k < n_target:
                zp = ping(zp, offsets(pat, k, h, prng), amp, d, b)
                k += 1
                for _ in range(interval):
                    ep, eu = shared(noise(nrng, zp.shape), zp, zu)
                    zp, zu = rk4(zp, MU, beta) + ep, rk4(zu, MU, beta) + eu
            # relax 16 substeps on copies (do not disturb the running sequence)
            rp, ru = zp.copy(), zu.copy()
            r2 = np.random.default_rng(900 + sd + n_target)
            for _ in range(16):
                ep, eu = shared(noise(r2, rp.shape), rp, ru)
                rp, ru = rk4(rp, MU, beta) + ep, rk4(ru, MU, beta) + eu
            res[n_target] = stats(shift_of(rp, ru, d, b), h)
        out[pat] = res
    return out


# ---------- Part B: navigating --------------------------------------------------------------
def part_b(sd, amp, every=2):
    _, d, b = g2.bank(sd)
    x, v = g2.trajectories(np.random.default_rng(300 + sd), E_B)
    out = {}
    for pat in PATTERNS:
        zp = np.sqrt(MU) * np.exp(1j * b[:, None] * (d @ x[0].T))
        zu = zp.copy()
        nrng = np.random.default_rng(70 + sd)
        prng = np.random.default_rng(80 + sd)
        k = 0
        hsum = np.zeros((E_B, 2))
        for t in range(STEPS):
            lin = MU + 1j * (b[:, None] * (d @ v[t].T) / (SUB * DT))
            if t % every == 0:
                sp = np.linalg.norm(v[t], axis=1, keepdims=True)
                h = v[t] / np.maximum(sp, 1e-12)
                hsum += h
                zp = ping(zp, offsets(pat, k, h, prng), amp, d, b)
                k += 1
            for _ in range(SUB):
                ep, eu = shared(noise(nrng, zp.shape), zp, zu)
                zp, zu = rk4(zp, lin, 0.0) + ep, rk4(zu, lin, 0.0) + eu
        sh = shift_of(zp, zu, d, b)
        hbar = hsum / np.maximum(np.linalg.norm(hsum, axis=1, keepdims=True), 1e-12)
        st = stats(sh, hbar)
        # the forward part could in principle be absorbed by a speed-gain recalibration:
        # how much of the shift lies along the total path travelled?
        path = x[-1] - x[0]
        g = (sh * path).sum() / max((path * path).sum(), 1e-12)
        resid = sh - g * path
        st["gain_absorbable_fraction"] = float(1 - (resid ** 2).sum() / max((sh ** 2).sum(), 1e-12))
        st["shift_after_gain_fix"] = float(np.linalg.norm(resid, axis=1).mean())
        st["sweeps"] = k
        out[pat] = st
    return out


def mean_sd(vals):
    return [float(np.mean(vals)), float(np.std(vals))]


def main():
    t0 = time.time()
    A = {}
    for beta in (0.0, 0.5):
        for amp in (0.1, 0.3):
            for interval in (8, 2):
                key = f"beta{beta:g}_a{amp:g}_int{interval}"
                per = [part_a(sd, beta, amp, interval) for sd in SEEDS]
                A[key] = {p: {str(n): {m: mean_sd([q[p][n][m] for q in per]) for m in per[0][p][n]}
                              for n in NS} for p in PATTERNS}
                print("A", key, {p: round(A[key][p]["16"]["mean_abs_shift"][0], 5) for p in PATTERNS},
                      f"{time.time() - t0:.0f}s", flush=True)
    B = {}
    for amp in (0.1, 0.3):
        per = [part_b(sd, amp) for sd in SEEDS]
        B[f"a{amp:g}"] = {p: {m: mean_sd([q[p][m] for q in per]) for m in per[0][p]} for p in PATTERNS}
        print("B", amp, {p: round(B[f"a{amp:g}"][p]["mean_abs_shift"][0], 5) for p in PATTERNS},
              f"{time.time() - t0:.0f}s", flush=True)

    def ratio(tab, p):
        return tab["same"]["mean_abs_shift"][0] / max(tab[p]["mean_abs_shift"][0], 1e-12)

    a0 = {p: A["beta0_a0.3_int8"][p]["16"] for p in PATTERNS}
    a5 = {p: A["beta0.5_a0.3_int8"][p]["16"] for p in PATTERNS}
    b3 = B["a0.3"]
    V = {"H1_left_right_halves_cost": bool(ratio(a0, "left_right") >= 2 and ratio(b3, "left_right") >= 2),
         "H1_ratios_same_over_lr": [ratio(a0, "left_right"), ratio(b3, "left_right")],
         "H2_antipodal_10x": bool(ratio(a0, "antipodal") >= 10 and ratio(b3, "antipodal") >= 10),
         "H2_ratios_same_over_antipodal": [ratio(a0, "antipodal"), ratio(b3, "antipodal")],
         "H3_shear_breaks_antipodal": bool(ratio(a5, "antipodal") < 10),
         "H3_ratio_beta0.5": ratio(a5, "antipodal")}
    out = {"setup": {"N": N, "sweep_length": S_LEN, "angle_deg": 30, "seeds": SEEDS, "E_A": E_A, "E_B": E_B,
                     "n_pings": NS, "patterns": PATTERNS, "partB_sweep_every_steps": 2},
           "partA": A, "partB": B, "verdicts": V, "seconds": round(time.time() - t0, 1)}
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/gate5_receipt.json", "w"), indent=2)
    print(json.dumps(V, indent=2), out["seconds"])


if __name__ == "__main__":
    main()
