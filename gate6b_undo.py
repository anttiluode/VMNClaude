"""Gate 6b - why the 'exact' conjugate compensator was no better than a sign flip.

Diagnostic run after Gate 6. The conjugate kick a*u^2*e^{-i phi} takes u from the bank's CURRENT
phase, which the first ping has already moved by f(alpha). Its relative phase is therefore
-(alpha - f(alpha)), not -alpha, which leaves a second-order residual, the same order as the
sign flip. Exact undo needs the PRE-QUERY phase u0: kick a*u*u0*e^{-i phi} sits at relative -alpha
and cancels to all orders (arg(1 + a e^{i a}) is odd).

Also: the same bank projected onto its limit cycle with noise off, to test whether the
pre-query undo is exact when the stored amplitudes are exactly sqrt(mu).

Same single bank and goal as one Gate 6 episode, noise on (sigma = 0.05) and off, listen window
16 and 64 substeps. Reports retained RMS per-unit phase disturbance against the unqueried twin.
"""
import json

import numpy as np

import gate2_path as g2
import gate4_ping as g4
import gate5_sweeps as g5
import gate6_restore as g6


def run(sd, Z, q, scheme, amp, window, rng):
    _, d, b = g2.bank(sd)
    phi = b[:, None] * (d @ q.T)
    zp, zu = Z.copy(), Z.copy()
    u0 = zp / np.abs(zp)
    kicks = [lambda z: amp * np.sqrt(g5.MU) * np.exp(1j * phi)]
    if scheme == "sign_flip":
        kicks.append(lambda z: -amp * np.sqrt(g5.MU) * np.exp(1j * phi))
    elif scheme == "conjugate_current":
        kicks.append(lambda z: amp * np.sqrt(g5.MU) * (z / np.abs(z)) ** 2 * np.exp(-1j * phi))
    elif scheme == "conjugate_prequery":
        kicks.append(lambda z: amp * np.sqrt(g5.MU) * (z / np.abs(z)) * u0 * np.exp(-1j * phi))
    for k in kicks:
        zp = zp + k(zp)
        for _ in range(window):
            zp, zu = g6.step(zp, zu, rng)
    for _ in range(g6.RELAX):
        zp, zu = g6.step(zp, zu, rng)
    return float(np.sqrt((np.angle(zp * np.conj(zu)) ** 2).mean()))


def main():
    sd = 10
    x, v = g2.trajectories(np.random.default_rng(310), 200)
    Z = g4.integrate(sd, "none", 0.0, x, v, np.random.default_rng(1))
    rng0 = np.random.default_rng(2)
    ang = rng0.uniform(0, 2 * np.pi, 200); rad = 0.35 * np.sqrt(rng0.uniform(0, 1, 200))
    q = x[-1] + np.stack([rad * np.cos(ang), rad * np.sin(ang)], 1)
    out = {}
    sigma0 = g5.SIGMA
    for sig in (sigma0, 0.0):
        g5.SIGMA = sig
        for window in (16, 64):
            for amp in (0.03, 0.1, 0.3):
                key = f"sigma{sig:g}_window{window}_a{amp:g}"
                out[key] = {s: run(sd, Z, q, s, amp, window, np.random.default_rng(0))
                            for s in ("single", "sign_flip", "conjugate_current", "conjugate_prequery")}
                print(key, {k: f"{v:.1e}" for k, v in out[key].items()}, flush=True)
    # same bank projected onto its limit cycle (|z| = sqrt(mu)), noise off: is the pre-query undo exact?
    g5.SIGMA = 0.0
    Zn = np.sqrt(g5.MU) * Z / np.abs(Z)
    out["end_amplitude_over_sqrt_mu_mean_sd"] = [float((np.abs(Z) / np.sqrt(g5.MU)).mean()), float((np.abs(Z) / np.sqrt(g5.MU)).std())]
    for amp in (0.03, 0.3):
        key = f"on_limit_cycle_sigma0_window64_a{amp:g}"
        out[key] = {s: run(sd, Zn, q, s, amp, 64, np.random.default_rng(0))
                    for s in ("single", "sign_flip", "conjugate_current", "conjugate_prequery")}
        print(key, {k: f"{v:.1e}" for k, v in out[key].items()}, flush=True)
    g5.SIGMA = sigma0
    json.dump(out, open("results/gate6b_receipt.json", "w"), indent=2)


if __name__ == "__main__":
    main()
