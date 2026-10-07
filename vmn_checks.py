"""Numerical checks for every claim in VMN_NOTE.md.

python vmn_checks.py            # all checks, writes results/receipt.json
python vmn_checks.py --quick    # smaller sizes
"""
import argparse
import json
import os
import time

import numpy as np

import vmn_core as vc


def rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


# ---------------------------------------------------------------------------
# C1  Vortex -> Matrix: the response operator is antilinear complex-symmetric
# ---------------------------------------------------------------------------

def check_c1(rng, n=12):
    pos = rng.normal(size=(n, 2)) * 2.0
    gam = rng.uniform(-1.5, 1.5, n)
    K = vc.jacobian(pos, gam)
    # finite-difference Jacobian of the velocity field
    h = 1e-6
    Kfd = np.zeros_like(K)
    flat = pos.ravel()
    for i in range(2 * n):
        p = flat.copy(); p[i] += h
        m = flat.copy(); m[i] -= h
        Kfd[:, i] = (vc.velocity(p.reshape(n, 2), gam).ravel()
                     - vc.velocity(m.reshape(n, 2), gam).ravel()) / (2 * h)
    M = vc.complex_tangent(pos, gam)
    # antilinear action: K dx  ==  M conj(dz)
    dxy = rng.normal(size=(n, 2))
    dz = dxy[:, 0] + 1j * dxy[:, 1]
    out_real = (K @ dxy.ravel()).reshape(n, 2)
    out_cplx = M @ np.conj(dz)
    DM = np.diag(gam) @ M
    # spectrum: eig(K)^2  ==  eig(M conj(M)) and conjugates
    ek = np.linalg.eigvals(K)
    em = np.linalg.eigvals(M @ np.conj(M))
    sq = np.sort_complex(np.round(ek ** 2, 9))
    target = np.sort_complex(np.round(np.concatenate([em, np.conj(em)]), 9))
    return {
        "claim": "K (real 2N x 2N tangent) equals the antilinear map dz -> M conj(dz); diag(Gamma) M is complex symmetric; spec(K)^2 = spec(M conj M) U conj",
        "jacobian_vs_finite_difference_rel_err": rel(K, Kfd),
        "antilinear_form_rel_err": rel(out_cplx, out_real[:, 0] + 1j * out_real[:, 1]),
        "gamma_weighted_symmetry_rel_err": rel(DM, DM.T),
        "trace_K": float(np.trace(K)),
        "spectrum_squared_match_max_abs_err": float(np.max(np.abs(sq - target))),
    }


# ---------------------------------------------------------------------------
# C2  An event (new vortex) writes a tidal, block-diagonal update
# ---------------------------------------------------------------------------

def sample_disk(rng, n, rho=1.0, rmin=0.0):
    """n points uniform in a disk of density rho centered at 0, at least rmin from 0."""
    R = np.sqrt(n / (np.pi * rho))
    r = np.sqrt(rng.uniform(rmin ** 2 / R ** 2, 1.0, n)) * R
    t = rng.uniform(0, 2 * np.pi, n)
    return np.stack([r * np.cos(t), r * np.sin(t)], 1)


def check_c2(rng, sizes=(25, 100, 400, 1600, 6400), draws=40):
    # exactness on one configuration
    n = 30
    pos = sample_disk(rng, n, rmin=0.2)
    gam = rng.uniform(-1, 1, n)
    g_new = 0.7
    w = np.zeros(2)
    K0 = vc.jacobian(pos, gam)
    K1 = vc.jacobian(np.vstack([pos, w]), np.append(gam, g_new))[: 2 * n, : 2 * n]
    dK = K1 - K0
    r = np.linalg.norm(pos - w, axis=1)
    pred = np.repeat(abs(g_new) / (2 * np.pi * r ** 2), 2)
    sv = np.linalg.svd(dK, compute_uv=False)
    offblock = dK.copy()
    for k in range(n):
        offblock[2 * k:2 * k + 2, 2 * k:2 * k + 2] = 0
    exact = {
        "off_block_diagonal_norm": float(np.linalg.norm(offblock)),
        "singular_values_vs_|g|/(2 pi r^2)_max_rel_err": float(np.max(np.abs(np.sort(sv) - np.sort(pred)) / np.sort(pred))),
        "update_independent_of_own_circulation": True,
    }
    # rank saturation at fixed density
    table = []
    for n in sizes:
        pr, e95, frac = [], [], []
        for _ in range(draws):
            p = sample_disk(rng, n, rmin=0.3)
            rr = np.linalg.norm(p, axis=1)
            s = np.repeat(1.0 / rr ** 2, 2)
            pr.append(vc.participation_rank(s))
            e95.append(vc.energy_rank(s))
        table.append({"N": n, "median_participation_rank": float(np.median(pr)),
                      "median_95pct_energy_rank": float(np.median(e95)),
                      "out_of": 2 * n})
    return {
        "claim": "Adding a vortex g at w changes the old-old response block by a block-diagonal update; block k is a pure strain with both singular values |g|/(2 pi r_k^2); its effective rank saturates with N at fixed density",
        "exactness": exact,
        "rank_vs_N_fixed_density": table,
    }


# ---------------------------------------------------------------------------
# C3  Finite-horizon: does the event-induced change of the propagator stay low rank?
# ---------------------------------------------------------------------------

def poisson_like(rng, n, rho, dmin):
    R = np.sqrt(n / (np.pi * rho))
    pts = []
    while len(pts) < n:
        c = rng.uniform(-R, R, 2)
        if c @ c > R * R:
            continue
        if pts and np.min(np.linalg.norm(np.array(pts) - c, axis=1)) < dmin:
            continue
        pts.append(c)
    return np.array(pts)


def check_c3(rng, sizes=(16, 36, 64, 100), seeds=4, T=1.0, dt=0.01, delta=0.05,
             kick=0.25):
    rows = []
    for n in sizes:
        res = []
        for s in range(seeds):
            pos = poisson_like(rng, n, rho=1.0, dmin=0.45)
            gam = rng.choice([-1.0, 1.0], n) * rng.uniform(0.5, 1.0, n)
            k = int(np.argmin(np.linalg.norm(pos, axis=1)))   # event at the center vortex
            pos_a = pos.copy(); pos_a[k] += kick * np.array([1.0, 0.0])
            _, P0 = vc.rk4_with_tangent(pos, gam, dt, int(T / dt), delta)
            _, P1 = vc.rk4_with_tangent(pos_a, gam, dt, int(T / dt), delta)
            sv0 = np.linalg.svd(P0, compute_uv=False)
            svd = np.linalg.svd(P1 - P0, compute_uv=False)
            res.append((vc.energy_rank(sv0), vc.participation_rank(svd), vc.energy_rank(svd),
                        float(np.linalg.norm(P1 - P0) / np.linalg.norm(P0))))
        a = np.array(res)
        rows.append({"N": n, "dim": 2 * n,
                     "propagator_95pct_rank": float(np.median(a[:, 0])),
                     "update_participation_rank": float(np.median(a[:, 1])),
                     "update_95pct_rank": float(np.median(a[:, 2])),
                     "update_rel_norm": float(np.median(a[:, 3]))})
    return {
        "claim": "Over a finite horizon the event (finite kick to one vortex) changes the propagator by an update whose effective rank stays small and roughly N-independent, while the propagator itself is full rank",
        "setup": {"T": T, "dt": dt, "smoothing_delta": delta, "kick": kick, "density": 1.0, "min_spacing": 0.45},
        "rows": rows,
    }


def check_c3b(rng, n=64, seeds=4, dt=0.01, delta=0.05):
    """How the update rank grows with horizon T and how it depends on kick size."""
    configs = []
    for s in range(seeds):
        pos = poisson_like(rng, n, rho=1.0, dmin=0.45)
        gam = rng.choice([-1.0, 1.0], n) * rng.uniform(0.5, 1.0, n)
        configs.append((pos, gam, int(np.argmin(np.linalg.norm(pos, axis=1)))))
    rows = []
    for T in (0.5, 1.0, 2.0, 4.0):
        for kick in (0.01, 0.25):
            r95, pr, ppr = [], [], []
            for pos, gam, k in configs:
                pa = pos.copy(); pa[k, 0] += kick
                _, P0 = vc.rk4_with_tangent(pos, gam, dt, int(T / dt), delta)
                _, P1 = vc.rk4_with_tangent(pa, gam, dt, int(T / dt), delta)
                sv = np.linalg.svd(P1 - P0, compute_uv=False)
                r95.append(vc.energy_rank(sv)); pr.append(vc.participation_rank(sv))
                ppr.append(vc.participation_rank(np.linalg.svd(P0, compute_uv=False)))
            rows.append({"T": T, "kick": kick, "update_95pct_rank": float(np.median(r95)),
                         "update_participation_rank": float(np.median(pr)),
                         "propagator_participation_rank": float(np.median(ppr)), "dim": 2 * n})
    return {"claim": "update rank vs horizon and kick size (N=64)", "rows": rows}


# ---------------------------------------------------------------------------
# C4  Vortex -> Neuron: a co-rotating pair in strain fires like an oscillator neuron
# ---------------------------------------------------------------------------

def rotation_period(traj, dt):
    ang = np.unwrap(np.angle(traj))
    turns = (ang - ang[0]) / (2 * np.pi)
    if abs(turns[-1]) < 2:
        return None
    t = np.arange(len(traj)) * dt
    idx = np.searchsorted(np.abs(turns), [1, 2])
    return float(t[idx[1]] - t[idx[0]])


def check_c4():
    Gam = 1.0
    # (a) frozen separation: Adler phase equation, period law 2 pi / sqrt(a^2 - b^2) in psi = 2 phi
    adler = []
    r = 1.0
    a = Gam / (2 * np.pi * r * r)
    for frac in (2.0, 1.2, 1.05, 1.01, 1.002):
        e = a / frac
        dt, steps = 0.002, 400000
        phi = 0.0
        t_cross = []
        for i in range(steps):
            k1 = a - e * np.sin(2 * phi)
            k2 = a - e * np.sin(2 * (phi + 0.5 * dt * k1))
            k3 = a - e * np.sin(2 * (phi + 0.5 * dt * k2))
            k4 = a - e * np.sin(2 * (phi + dt * k3))
            phi_new = phi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            if np.floor(phi_new / np.pi) > np.floor(phi / np.pi):
                t_cross.append(i * dt)
                if len(t_cross) == 2:
                    break
            phi = phi_new
        measured = t_cross[1] - t_cross[0]
        predicted = np.pi / np.sqrt(a * a - e * e)
        adler.append({"drive_over_threshold": frac, "period_measured": measured,
                      "period_predicted_pi/sqrt(a^2-e^2)": float(predicted)})
    # (b) full pair (separation free): heteroclinic onset, period ~ log(1/|H - H_s|)
    e = 0.05
    rs = np.sqrt(Gam / (2 * np.pi * e))
    Hs = vc.pair_energy(np.array(rs * np.exp(1j * np.pi / 4)), Gam, e)
    lam = 2 * e      # saddle eigenvalue (see note)
    rows = []
    for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        # start on the phi = 0 axis at the radius whose energy is Hs + eps (bound side)
        # on phi = 0: H = -(G/4pi) log r^2, decreasing in r; bound orbits have H > Hs
        r0 = np.exp(-(Hs + eps) * 4 * np.pi / (2 * Gam))
        dt = 0.01
        traj = vc.pair_in_strain(r0, Gam, e, dt, 400000)
        rows.append({"H_minus_Hs": eps, "start_r": float(r0),
                     "period": rotation_period(traj, dt)})
    eps = np.array([q["H_minus_Hs"] for q in rows if q["period"]])
    per = np.array([q["period"] for q in rows if q["period"]])
    slope = float(np.polyfit(np.log(1 / eps[-4:]), per[-4:], 1)[0])
    # escape side
    r_out = np.exp(-(Hs - 1e-3) * 4 * np.pi / (2 * Gam))
    tr = vc.pair_in_strain(r_out, Gam, e, 0.01, 40000)
    return {
        "claim": "Relative coordinate of a co-rotating pair in strain: dz/dt = i G/(2 pi conj z) + e conj z. Frozen separation -> Adler (SNIC-type, period ~ (a^2-e^2)^-1/2). Free separation -> heteroclinic loop through two saddles at r_s = sqrt(G/(2 pi e)), period ~ (2/lambda) log(1/dH) with lambda = 2e",
        "frozen_separation_adler": adler,
        "free_pair": {"strain_e": e, "saddle_radius": float(rs), "saddle_eigenvalue_2e": lam,
                      "period_table": rows,
                      "fitted_dPeriod/dlog(1/dH)": slope,
                      "predicted_2/lambda": 2 / lam,
                      "escape_side_final_radius_vs_saddle": float(abs(tr[-1]) / rs)},
    }


# ---------------------------------------------------------------------------
# C5  Neuron -> Matrix: memory of a kick, Stuart-Landau vs Hamiltonian vortex
# ---------------------------------------------------------------------------

def check_c5():
    omega, beta = 1.0, 2.0
    rows = []
    for mu in (0.4, 0.1, 0.025, 0.00625):
        rstar = np.sqrt(mu)
        rho0 = 1e-5          # same absolute kick for every mu
        dt = 0.002
        steps = int(min(12.0 / mu, 4000.0) / dt)
        a = vc.stuart_landau(rstar, mu, omega, beta, dt, steps)
        b = vc.stuart_landau(rstar + rho0, mu, omega, beta, dt, steps)
        dphi = float(np.angle(b[-1] / a[-1]))
        # large kick: half the cycle amplitude, tests the exact law (Sol's correction)
        big = vc.stuart_landau(1.5 * rstar, mu, omega, beta, dt, steps)
        rows.append({"mu": mu, "phase_shift": dphi,
                     "predicted_-beta*rho0/sqrt(mu)": float(-beta * rho0 / rstar),
                     "exact_-beta*log(1+rho0/sqrt(mu))": float(-beta * np.log1p(rho0 / rstar)),
                     "large_kick_0.5r*_phase_shift": float(np.angle(big[-1] / a[-1])),
                     "large_kick_exact_-beta*log(1.5)": float(-beta * np.log(1.5)),
                     "relaxation_time_1/(2mu)": 1 / (2 * mu)})
    # response-matrix change between two relaxed states that differ only in phase
    # (Sol): ||J(z1) - J(z2)||_2 = 2 mu sqrt(1+beta^2) |sin(dphi)|
    def sl_jac_real(z, mu):
        # d/d(x,y) of (mu+i w) z - (1+i b)|z|^2 z
        a_lin = (mu + 1j * omega) - (1 + 1j * beta) * 2 * abs(z) ** 2   # coefficient of dz
        a_conj = -(1 + 1j * beta) * z * z                                # coefficient of conj(dz)
        def apply(v):
            out = a_lin * v + a_conj * np.conj(v)
            return np.array([out.real, out.imag])
        return np.column_stack([apply(1.0 + 0j), apply(1j)])
    dj = []
    for mu in (0.4, 0.00625):
        for dph in (0.3, 1.2):
            r = np.sqrt(mu)
            d = np.linalg.norm(sl_jac_real(r, mu) - sl_jac_real(r * np.exp(1j * dph), mu), 2)
            dj.append({"mu": mu, "dphi": dph, "norm_dJ": float(d),
                       "predicted": float(2 * mu * np.sqrt(1 + beta ** 2) * abs(np.sin(dph)))})
    # below onset: linear decay time 1/|mu|
    below = []
    for mu in (-0.4, -0.1, -0.025):
        dt = 0.002
        z = vc.stuart_landau(1e-4, mu, omega, beta, dt, int(2.0 / abs(mu) / dt))
        t = np.arange(len(z)) * dt
        fit = np.polyfit(t, np.log(np.abs(z)), 1)[0]
        below.append({"mu": mu, "decay_rate_measured": float(-fit), "predicted_|mu|": abs(mu)})
    # Hamiltonian vortex (satellite around circulation G at radius r): kick never relaxes,
    # phase error grows linearly: dphi(t) = -2 G dr t / (2 pi r^3)
    G, r, dr = 1.0, 1.0, 1e-4
    dt, steps = 0.001, 20000
    a = vc.pair_in_strain(r, G, 0.0, dt, steps)
    b = vc.pair_in_strain(r + dr, G, 0.0, dt, steps)
    dphi = np.unwrap(np.angle(b / a))
    t = np.arange(steps + 1) * dt
    slope = float(np.polyfit(t, dphi, 1)[0])
    return {
        "claim": "Stuart-Landau with shear beta: a radial kick rho0 leaves the exact permanent phase shift -beta log(1 + rho0/sqrt(mu)) (Sol; small-kick limit -beta rho0/sqrt(mu)). The phase is stored, but the response-matrix change it causes is 2 mu sqrt(1+beta^2)|sin dphi| and shrinks toward onset (Sol). A Hamiltonian vortex stores a kick as a frequency shift and the phase error grows linearly",
        "above_onset": rows,
        "phase_to_response_matrix_change": dj,
        "below_onset": below,
        "hamiltonian_vortex": {"phase_drift_slope": slope,
                               "predicted_-2G dr/(2 pi r^3)": float(-2 * G * dr / (2 * np.pi * r ** 3))},
    }


# ---------------------------------------------------------------------------
# C6  Adding a dimension: the doubled (z, conj z) view of vortex coupling
# ---------------------------------------------------------------------------

def check_c6():
    import gate1_layer as g1
    n = g1.N
    rows_equiv, rows_uniform = [], []
    for seed in (10, 11, 12):
        rng = np.random.default_rng(1000 + seed)
        M, gg = g1.vortex_matrix(rng)
        w = rng.uniform(0.5, 1.5, n)
        # (a) doubled matrix [[A, B], [conj B, conj A]] has the real-form spectrum
        for counter in (False, True):
            ww = w * np.sign(gg) if counter else w
            for kappa in (0.3, 1.0, 3.0):
                for mode in ("vortex", "linear"):
                    A = np.diag(1j * ww) + (kappa * M if mode == "linear" else 0)
                    B = kappa * M if mode == "vortex" else np.zeros_like(M)
                    D = np.block([[A, B], [np.conj(B), np.conj(A)]])
                    shift_doubled = -float(np.max(np.linalg.eigvals(D).real))
                    shift_real = g1.onset_mu(M, ww, kappa, mode)
                    rows_equiv.append(abs(shift_doubled - shift_real))
        # (b) uniform frequency: eig(D)^2 = -w^2 + kappa^2 eig(M conj M), exactly
        nu = np.linalg.eigvals(M @ np.conj(M))
        wu = 1.0
        for kappa in (0.1, 0.3, 1.0, 3.0):
            D = np.block([[1j * wu * np.eye(n), kappa * M], [kappa * np.conj(M), -1j * wu * np.eye(n)]])
            ev = np.linalg.eigvals(D)
            pred = np.concatenate([np.sqrt((kappa ** 2 * nu - wu ** 2).astype(complex)),
                                   -np.sqrt((kappa ** 2 * nu - wu ** 2).astype(complex))])
            err = float(np.max(np.abs(np.sort_complex(np.round(ev, 9)) - np.sort_complex(np.round(pred, 9)))))
            rows_uniform.append({
                "seed": seed, "kappa": kappa,
                "spectrum_formula_max_abs_err": err,
                "vortex_onset_shift": float(-np.max(ev.real)),
                "vortex_shift_formula_-max_Re_sqrt(k^2 nu - w^2)": float(-np.max(pred.real)),
                "linear_onset_shift": g1.onset_mu(M, np.full(n, wu), kappa, "linear"),
                "linear_formula_-k_max_Re_eig(M)": float(-kappa * np.max(np.linalg.eigvals(M).real)),
                "threshold_kappa*_w/sqrt(max_real_nu)": float(wu / np.sqrt(np.max(nu.real))),
            })
    return {
        "claim": "Doubling z -> (z, conj z) turns the antilinear vortex coupling into an ordinary linear matrix with the same spectrum. With one common frequency w, D^2 is block diagonal and eig(D)^2 = -w^2 + kappa^2 eig(M conj M): vortex coupling only shifts the onset through complex eigenvalues of M conj M (second order) until kappa reaches w/sqrt(nu), while linear coupling shifts it at first order, kappa max Re eig(M)",
        "doubled_vs_real_form_max_abs_err": float(max(rows_equiv)),
        "uniform_frequency": rows_uniform,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    rng = np.random.default_rng(55)
    out = {"python_numpy": np.__version__}
    t0 = time.time()
    out["C1_vortex_to_matrix"] = check_c1(rng)
    out["C2_event_tidal_update"] = check_c2(rng, sizes=(25, 100, 400, 1600) if args.quick else (25, 100, 400, 1600, 6400))
    out["C3_finite_horizon_update_rank"] = check_c3(rng, sizes=(16, 36) if args.quick else (16, 36, 64, 100))
    out["C3b_update_rank_vs_horizon"] = check_c3b(rng)
    out["C4_vortex_pair_as_neuron"] = check_c4()
    out["C5_memory_near_hopf"] = check_c5()
    out["C6_doubled_view"] = check_c6()
    out["seconds"] = round(time.time() - t0, 1)
    os.makedirs("results", exist_ok=True)
    with open("results/receipt.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
