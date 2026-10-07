"""Core objects for Vortex - Matrix - Neuron.

Point vortices in the plane, their exact tangent (response) operator in real
and complex form, and the two small oscillator models used as "neurons".
Everything is plain NumPy so each identity in VMN_NOTE.md can be checked.
"""
import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Point vortices
# ---------------------------------------------------------------------------

def velocity(pos, gam, delta=0.0):
    """Velocity of every vortex. pos: (N,2), gam: (N,). Krasny smoothing delta."""
    d = pos[:, None, :] - pos[None, :, :]            # d[k,j] = x_k - x_j
    s = (d ** 2).sum(-1) + delta ** 2
    np.fill_diagonal(s, np.inf)
    w = gam[None, :] / (TWO_PI * s)
    u = -(w * d[..., 1]).sum(1)
    v = (w * d[..., 0]).sum(1)
    return np.stack([u, v], 1)


def pair_block(d, g, delta=0.0):
    """d(velocity of k)/d(x_k - x_j) for one source of strength g at offset d.

    Returns the 2x2 block. For delta=0 it is symmetric and traceless:
    (g / 2 pi r^4) [[2xy, y^2-x^2], [y^2-x^2, -2xy]] -> singular values g/(2 pi r^2).
    """
    dx, dy = d
    s = dx * dx + dy * dy + delta ** 2
    c = g / (TWO_PI * s * s)
    return c * np.array([[2 * dx * dy, 2 * dy * dy - s],
                         [s - 2 * dx * dx, -2 * dx * dy]])


def jacobian(pos, gam, delta=0.0):
    """Exact real tangent operator K = d(velocity)/d(positions), shape (2N,2N)."""
    n = len(gam)
    d = pos[:, None, :] - pos[None, :, :]
    s = (d ** 2).sum(-1) + delta ** 2
    np.fill_diagonal(s, 1.0)                           # self term masked below
    c = gam[None, :] / (TWO_PI * s * s)                # (N,N)
    np.fill_diagonal(c, 0.0)
    dx, dy = d[..., 0], d[..., 1]
    b00 = c * 2 * dx * dy
    b01 = c * (2 * dy * dy - s)
    b10 = c * (s - 2 * dx * dx)
    b11 = -c * 2 * dx * dy
    for b in (b00, b01, b10, b11):
        np.fill_diagonal(b, 0.0)
    K = np.zeros((2 * n, 2 * n))
    # off-diagonal: d vel_k / d x_j = -B_kj
    K[0::2, 0::2] = -b00
    K[0::2, 1::2] = -b01
    K[1::2, 0::2] = -b10
    K[1::2, 1::2] = -b11
    # diagonal blocks: d vel_k / d x_k = sum_j B_kj
    idx = np.arange(n)
    K[2 * idx, 2 * idx] = b00.sum(1)
    K[2 * idx, 2 * idx + 1] = b01.sum(1)
    K[2 * idx + 1, 2 * idx] = b10.sum(1)
    K[2 * idx + 1, 2 * idx + 1] = b11.sum(1)
    return K


def complex_tangent(pos, gam):
    """M such that  d(dz)/dt = M conj(dz)  (unsmoothed point vortices).

    M_kj = (i/2pi) G_j / conj(z_k - z_j)^2,  M_kk = -sum_j M_kj.
    diag(G) M is complex symmetric, so the response is an antilinear
    complex-symmetric map (Takagi structure).
    """
    z = pos[:, 0] + 1j * pos[:, 1]
    dz = np.conj(z[:, None] - z[None, :])
    np.fill_diagonal(dz, 1.0)                           # self term overwritten below
    M = (1j / TWO_PI) * gam[None, :] / dz ** 2
    np.fill_diagonal(M, 0.0)
    M[np.diag_indices_from(M)] = -M.sum(1)
    return M


def rk4_with_tangent(pos, gam, dt, steps, delta):
    """Integrate vortices and the variational equation  dPhi/dt = K Phi."""
    n = len(gam)
    x = pos.copy()
    P = np.eye(2 * n)

    def f(x, P):
        return velocity(x, gam, delta), jacobian(x, gam, delta) @ P

    for _ in range(steps):
        k1x, k1P = f(x, P)
        k2x, k2P = f(x + 0.5 * dt * k1x, P + 0.5 * dt * k1P)
        k3x, k3P = f(x + 0.5 * dt * k2x, P + 0.5 * dt * k2P)
        k4x, k4P = f(x + dt * k3x, P + dt * k3P)
        x = x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        P = P + dt / 6 * (k1P + 2 * k2P + 2 * k3P + k4P)
    return x, P


# ---------------------------------------------------------------------------
# Rank measures
# ---------------------------------------------------------------------------

def participation_rank(sv):
    """(sum s^2)^2 / sum s^4 over singular values s."""
    e = np.asarray(sv) ** 2
    return float(e.sum() ** 2 / (e ** 2).sum())


def energy_rank(sv, frac=0.95):
    e = np.sort(np.asarray(sv) ** 2)[::-1]
    c = np.cumsum(e) / e.sum()
    return int(np.searchsorted(c, frac) + 1)


# ---------------------------------------------------------------------------
# Neurons
# ---------------------------------------------------------------------------

def pair_in_strain(z0, Gam, e, dt, steps):
    """Relative coordinate of two co-rotating vortices (total circulation Gam)
    in pure strain e:  dz/dt = i Gam / (2 pi conj z) + e conj z."""
    def f(z):
        return 1j * Gam / (TWO_PI * np.conj(z)) + e * np.conj(z)
    z = complex(z0)
    out = np.empty(steps + 1, complex)
    out[0] = z
    for i in range(steps):
        k1 = f(z); k2 = f(z + 0.5 * dt * k1)
        k3 = f(z + 0.5 * dt * k2); k4 = f(z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = z
    return out


def pair_energy(z, Gam, e):
    x, y = z.real, z.imag
    return -(Gam / (4 * np.pi)) * np.log(x * x + y * y) + e * x * y


def stuart_landau(z0, mu, omega, beta, dt, steps):
    """dz/dt = (mu + i omega) z - (1 + i beta) |z|^2 z."""
    def f(z):
        return (mu + 1j * omega) * z - (1 + 1j * beta) * abs(z) ** 2 * z
    z = complex(z0)
    out = np.empty(steps + 1, complex)
    out[0] = z
    for i in range(steps):
        k1 = f(z); k2 = f(z + 0.5 * dt * k1)
        k3 = f(z + 0.5 * dt * k2); k4 = f(z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = z
    return out
