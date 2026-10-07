# VMNClaude — Vortex · Matrix · Neuron

Claude (Opus 5.5) side of the VMN line. Sol works the same question in [VMN](https://github.com/anttiluode/VMN); several corrections here came from Sol's review and are credited where they appear.

**Question:** is there exact mathematics linking a vortex, a matrix and a neuron, with the nonlinearity doing the real work?

**Short answer, from the note and its checks:** yes, through one term. The nonlinearity that matters is **shear**, frequency that depends on amplitude. It appears as a vortex's self-rotation, as the β of an oscillator neuron, and as the amplitude-dependent wave speed in [Kompressori](https://github.com/anttiluode/Kompressori). Without a nonlinearity, no event can change how a system responds to the next input (§0).

Full derivations: **[VMN_NOTE.md](VMN_NOTE.md)**. First layer test: **[GATE1.md](GATE1.md)**. Every number is in a receipt under `results/`.

## Latest: Gate 1 — the layer, against kill conditions written in advance

![gate 1](results/gate1_summary.png)

40 oscillator neurons coupled through the vortex matrix, against the same matrix used as ordinary coupling, no coupling, and an echo state network of the same state size.

- **Vortex coupling is useless when all oscillators turn the same way.** Its conjugate link pairs each unit with a partner turning the opposite way, so it averages out: it moves the onset by 0.002 where ordinary coupling moves it by 0.126.
- **It wakes up when they turn both ways,** as vortices of both signs do. Then it beats ordinary coupling on both tasks: memory capacity 31.0 vs 27.5, NARMA10 error 0.391 vs 0.533. This was predicted before the run, but as a follow-up after seeing the first result, not pre-registered.
- **An echo state network still wins** (memory capacity 32–37, NARMA10 0.334). As a drop-in reservoir, VMN fails its own kill condition.
- **The best operating point is on the stable side, close to onset, and consistency breaks just above onset,** as §4 predicted and the cylinder-wake study reported.
- Shear (β) did not help on these tasks.

## The mathematics

![summary](results/vmn_summary.png)

| | Link | Result | Status |
|---|---|---|---|
| §1 | Vortex → Matrix | Vortex response is an antilinear, circulation-weighted complex-symmetric map; the frozen Jacobian's spectrum² = eigenvalues of $`M\bar M`$ | Exact, verified to 1e-16 |
| §2 | Event → low-rank update | A new vortex changes each old vortex's response by a tidal strain $`\lvert g\rvert/2\pi r^2`$; exact rank is full (2N), effective rank ≈ 5.7 from N = 25 to 6400 | Exact one-step; finite horizon measured |
| §2 | Finite horizon | Update rank 4–7 while propagator rank grows with N (T = 1) | Measured; meaningless at long T, where everything collapses |
| §3 | Vortex → Neuron | Co-rotating pair in strain: frozen = Adler/SNIC neuron, free = heteroclinic (log) onset | Verified; slope 20.03 vs 20.0 |
| §4 | Neuron → Matrix | Stuart–Landau with shear keeps exactly $`-\beta\log(1+\rho_0/\sqrt\mu)`$ of a kick (Sol's sharpening); the phase is stored but the response change it leaves, $`2\mu\sqrt{1+\beta^2}\lvert\sin\Delta\varphi\rvert`$, shrinks toward onset (Sol) | Verified to ~1e-10. Optimum just *below* onset, where memory still fades: confirmed in Gate 1 |

## What this does not show

- Nothing here is a new physical law. The ingredients are classical; the contribution is putting them side by side with checks.
- The layer has only been tested as an untrained reservoir, on two standard benchmarks, at N = 40. Trained versions, moving geometry and rhythm-structured tasks are untested.
- The finite-horizon low-rank result holds only for horizons shorter than the chaotic stretching time.

## Run

```bash
pip install numpy matplotlib
python vmn_checks.py         # ~70 s
python make_figure.py
python gate1_layer.py        # ~13 min
python gate1b_followups.py   # ~8 min
python make_gate1_figure.py
```

Files: `vmn_core.py` (vortex dynamics, exact Jacobians, oscillators) · `vmn_checks.py` (maths checks) · `gate1_layer.py`, `gate1b_followups.py` (layer test) · `make_figure.py`, `make_gate1_figure.py` · `results/`.
