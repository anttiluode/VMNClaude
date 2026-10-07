# VMNClaude — Vortex · Matrix · Neuron

Claude (Opus 5.5) side of the VMN line. Sol works the same question in [VMN](https://github.com/anttiluode/VMN); several corrections here came from Sol's review and are credited where they appear.

**Question:** is there exact mathematics linking a vortex, a matrix and a neuron, with the nonlinearity doing the real work?

**Short answer, from the note and its checks:** yes, through one term. The nonlinearity that matters is **shear**, frequency that depends on amplitude. It appears as a vortex's self-rotation, as the β of an oscillator neuron, and as the amplitude-dependent wave speed in [Kompressori](https://github.com/anttiluode/Kompressori). Without a nonlinearity, no event can change how a system responds to the next input (§0).

Full derivations: **[VMN_NOTE.md](VMN_NOTE.md)**. Layer tests: **[GATE1.md](GATE1.md)** (reservoir benchmarks) and **[GATE2.md](GATE2.md)** (path integration). Every number is in a receipt under `results/`.

## The entorhinal link

The entorhinal cortex's best-known models are built from the same pieces as VMN: oscillators whose frequency depends on a variable, with phase as memory.

- **Grid cells as velocity-controlled oscillators.** In the oscillatory interference model, each oscillator's frequency rises above a baseline in proportion to running speed along a preferred direction; its phase against the baseline integrates distance travelled, and oscillators 60° apart interfere into a hexagonal grid ([Jeewajee et al. 2008](https://ucl.ac.uk/icn/sites/icn/files/jee08.pdf)). That is the same mathematics as §4's Hamiltonian vortex, where a kick changes the frequency and the phase becomes an integrator. The model's known weakness, phase drift, is the price §4 names.
- **Position as rotation.** [Gao, Xie, Wei, Zhu & Wu (NeurIPS 2021)](https://arxiv.org/abs/2006.10259) model grid cells as a vector rotated by a matrix depending on displacement, $`v(x+\Delta x)=M(\Delta x)\,v(x)`$, and get hexagons from a distance-preserving condition plus three plane waves at 120°. That is the real-matrix form of a complex oscillator bank, the representation §5b works in.
- **A frequency gradient sets the scales.** Entorhinal stellate cells' resonance frequency varies along the dorsal–ventral axis, as grid spacing does; knocking out the HCN1 channel flattens the frequency gradient ([Giocomo & Hasselmo 2009](https://www.bu.edu/hasselmo/GiocomoHasselmoJNeur2009.pdf)). In VMN terms, the spread of unit frequencies sets which spatial scales a bank can represent.
- **Hypothesis from Gate 1b (untested in the brain):** path integration contains both rotation senses, because moving with or against a preferred direction turns the relative phase opposite ways. That is exactly where conjugate (vortex) coupling becomes resonant. I found nothing on whether the brain uses anything like it.
- **The hexagon–vortex-lattice look-alike.** Grid hexagons resemble the hexagonal vortex lattices of superconductors and rotating superfluids, and the analogy has been drawn ([Chalyi, Ukr. J. Phys. 2022](https://ujp.bitp.kiev.ua/index.php/uik/en/article/view/2022537)). But hexagons appear whenever repelling objects pack in 2D, so the shared shape is not evidence of a shared mechanism.
- **AI cousin:** rotary position embeddings (RoPE) in transformers encode position by rotating vectors by angles proportional to position, which is the Gao et al. picture. The grid-cell resemblance has, I believe, been pointed out by others; not checked.

## Latest: Gate 3 — memory vs control, and event-gated coupling

![gate 3](results/gate3_summary.png)

Full write-up: **[GATE3.md](GATE3.md)**. Gate 1 (coupling helps) and Gate 2 (coupling hurts) are two ends of one trade-off, fixed by symmetry. An uncoupled bank's phases are perfect integrators, but the response operator cannot see what they store. Coupling lets the operator see the stored state (control) and makes the phases drift. Kill conditions were written before running, and all three pass.

- **The trade-off is real.** Control is exactly 0 without coupling and grows linearly with κ. No coupling strength lowers path-integration error; past κ ≈ 0.01, error climbs to chance by κ = 0.3.
- **Global phase:** rotating all phases together leaves the operator's singular values untouched under linear coupling (4×10⁻¹⁶) and changes them by 3.8% under vortex coupling.
- **Event gating wins:** coupling switched on only 10% of the time gives **7–9× more control** than steady coupling at the same error. That is close to the expected ~10×, so this confirms the scaling rather than surprising.
- **Not pre-registered:** vortex coupling beats linear on *both* axes at every κ ≥ 0.01, with more control and less error.
- **Not shown:** that the control is *useful*. That needs a task where the answer must depend on stored position, which is the next gate.

Design rule: keep the integrator symmetric, and break the symmetry briefly, only when a state-dependent response is needed.

## Gate 2 — path integration, the entorhinal setting

![gate 2](results/gate2_summary.png)

Full write-up: **[GATE2.md](GATE2.md)**. A bank of 40 velocity-controlled oscillators gets velocity only, after a start landmark, with noise on every unit. Kill conditions were written before running.

- **The entorhinal mechanism works:** end-of-episode position error 0.18, against 0.36 for an echo state network (0.38 is guessing the centre). The credit goes to velocity-controlled frequency, not to coupling: the *uncoupled* bank is best.
- **Vortex coupling adds nothing here** (K1 fires). Tuning drove all coupling to its minimum. Path integration needs each phase to stay the pure integral of its own input, so any coupling only adds error. Gate 1b's rotation-sense effect helps where mixing is the job, not here.
- **Shear hurts an integrator:** β = 2 nearly doubles the error (0.18 → 0.33), because shear turns amplitude noise into phase drift (§4). Integrating oscillators should be isochronous.
- A demodulation bug that made all β = 2 runs look like chance was caught and fixed before this write-up; GATE2.md says how.

## Gate 1 — the layer, against kill conditions written in advance

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
| §5b | Adding a dimension | Doubling z → (z, z̄) turns vortex coupling into an ordinary matrix; with one frequency $`\lambda^2=\kappa^2\nu-\omega^2`$, which is why vortex coupling cancelled in Gate 1 (second order in κ below a threshold, first order for ordinary coupling) | Exact, verified to 1e-14 |

## What this does not show

- Nothing here is a new physical law. The ingredients are classical; the contribution is putting them side by side with checks.
- The layer has only been tested as an untrained reservoir at N = 40: two standard benchmarks and path integration. Trained versions, moving geometry, direction-structured coupling and rhythm tasks are untested.
- The finite-horizon low-rank result holds only for horizons shorter than the chaotic stretching time.

## Run

```bash
pip install numpy matplotlib
python vmn_checks.py         # ~70 s
python make_figure.py
python gate1_layer.py        # ~13 min
python gate1b_followups.py   # ~8 min
python make_gate1_figure.py
python gate2_path.py         # ~5 min
python make_gate2_figure.py
```

Files: `vmn_core.py` (vortex dynamics, exact Jacobians, oscillators) · `vmn_checks.py` (maths checks) · `gate1_layer.py`, `gate1b_followups.py` (reservoir test) · `gate2_path.py` (path integration) · `make_figure.py`, `make_gate1_figure.py`, `make_gate2_figure.py` · `results/`.
