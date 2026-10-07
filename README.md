# VMNClaude — Vortex · Matrix · Neuron

Claude (Opus 5.5) side of the VMN line. Sol works the same question in its own repo.

**Question:** is there exact mathematics linking a vortex, a matrix and a neuron, with the nonlinearity doing the real work?

**Short answer, from the note and its checks:** yes, through one term. The nonlinearity that matters is **shear**, frequency that depends on amplitude. It appears as a vortex's self-rotation, as the β of an oscillator neuron, and as the amplitude-dependent wave speed in [Kompressori](https://github.com/anttiluode/Kompressori). Without a nonlinearity, no event can change how a system responds to the next input (§0).

Full derivations: **[VMN_NOTE.md](VMN_NOTE.md)**. Every number below is in `results/receipt.json`.

![summary](results/vmn_summary.png)

| | Link | Result | Status |
|---|---|---|---|
| §1 | Vortex → Matrix | Vortex response is an antilinear, circulation-weighted complex-symmetric map; frequencies² = eigenvalues of $`M\bar M`$ | Exact, verified to 1e-16 |
| §2 | Event → low-rank update | A new vortex changes each old vortex's response by a tidal strain $`\lvert g\rvert/2\pi r^2`$; effective rank ≈ 5.7 from N = 25 to 6400 | Exact one-step; finite horizon measured |
| §2 | Finite horizon | Update rank 4–7 while propagator rank grows with N (T = 1) | Measured; meaningless at long T, where everything collapses |
| §3 | Vortex → Neuron | Co-rotating pair in strain: frozen = Adler/SNIC neuron, free = heteroclinic (log) onset | Verified; slope 20.03 vs 20.0 |
| §4 | Neuron → Matrix | Stuart–Landau with shear keeps $`-\beta\rho_0/\sqrt\mu`$ of a kick; memory diverges at the Hopf point | Verified; as an explanation of the Kármán computing peak it is a hypothesis |

## What this does not show

- Nothing here is a new physical law. The ingredients are classical; the contribution is putting them side by side with checks.
- No AI layer has been built or tested yet. §6 of the note proposes one with kill conditions fixed in advance.
- The finite-horizon low-rank result holds only for horizons shorter than the chaotic stretching time.

## Run

```bash
pip install numpy matplotlib
python vmn_checks.py      # ~70 s
python make_figure.py
```

Files: `vmn_core.py` (vortex dynamics, exact Jacobians, oscillators) · `vmn_checks.py` (all checks) · `make_figure.py` · `results/`.
