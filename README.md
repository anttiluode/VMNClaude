# VMNClaude — Vortex · Matrix · Neuron

Claude (Opus 5.5) side of the VMN line. Sol works the same question in [VMN](https://github.com/anttiluode/VMN); several corrections here came from Sol's review and are credited where they appear.

**Question:** is there exact mathematics linking a vortex, a matrix and a neuron, with the nonlinearity doing the real work?

**Short answer, from the note and its checks:** yes, through one term. The nonlinearity that matters is **shear**, frequency that depends on amplitude. It appears as a vortex's self-rotation, as the β of an oscillator neuron, and as the amplitude-dependent wave speed in [Kompressori](https://github.com/anttiluode/Kompressori). Without a nonlinearity, no event can change how a system responds to the next input (§0).

Full derivations: **[VMN_NOTE.md](VMN_NOTE.md)**. Layer tests: **[GATE1.md](GATE1.md)** (reservoir benchmarks), **[GATE2.md](GATE2.md)** (path integration), **[GATE3.md](GATE3.md)** (memory vs control) and **[GATE4.md](GATE4.md)** (ping-and-listen goal vectors). Every number is in a receipt under `results/`.

## The entorhinal link

The entorhinal cortex's best-known models are built from the same pieces as VMN: oscillators whose frequency depends on a variable, with phase as memory.

- **Grid cells as velocity-controlled oscillators.** In the oscillatory interference model, each oscillator's frequency rises above a baseline in proportion to running speed along a preferred direction; its phase against the baseline integrates distance travelled, and oscillators 60° apart interfere into a hexagonal grid ([Jeewajee et al. 2008](https://ucl.ac.uk/icn/sites/icn/files/jee08.pdf)). That is the same mathematics as §4's Hamiltonian vortex, where a kick changes the frequency and the phase becomes an integrator. The model's known weakness, phase drift, is the price §4 names.
- **Position as rotation.** [Gao, Xie, Wei, Zhu & Wu (NeurIPS 2021)](https://arxiv.org/abs/2006.10259) model grid cells as a vector rotated by a matrix depending on displacement, $`v(x+\Delta x)=M(\Delta x)\,v(x)`$, and get hexagons from a distance-preserving condition plus three plane waves at 120°. That is the real-matrix form of a complex oscillator bank, the representation §5b works in.
- **A frequency gradient sets the scales.** Entorhinal stellate cells' resonance frequency varies along the dorsal–ventral axis, as grid spacing does; knocking out the HCN1 channel flattens the frequency gradient ([Giocomo & Hasselmo 2009](https://www.bu.edu/hasselmo/GiocomoHasselmoJNeur2009.pdf)). In VMN terms, the spread of unit frequencies sets which spatial scales a bank can represent.
- **Hypothesis from Gate 1b (untested in the brain):** path integration contains both rotation senses, because moving with or against a preferred direction turns the relative phase opposite ways. That is exactly where conjugate (vortex) coupling becomes resonant. I found nothing on whether the brain uses anything like it.
- **The hexagon–vortex-lattice look-alike.** Grid hexagons resemble the hexagonal vortex lattices of superconductors and rotating superfluids, and the analogy has been drawn ([Chalyi, Ukr. J. Phys. 2022](https://ujp.bitp.kiev.ua/index.php/uik/en/article/view/2022537)). But hexagons appear whenever repelling objects pack in 2D, so the shared shape is not evidence of a shared mechanism.
- **AI cousin:** rotary position embeddings (RoPE) in transformers encode position by rotating vectors by angles proportional to position, which is the Gao et al. picture. The grid-cell resemblance has, I believe, been pointed out by others; not checked.

## Does anyone need this? — [NEEDS.md](NEEDS.md)

An honest market check after Gate 4: **almost nothing needs ping-queryable oscillator memory more than a digital register, and where it is needed, it already exists.** VCO-based ADCs use oscillator phase as an integrator; passive SAW tags are read by ping; frequency-multiplexed resonator arrays are read on one wire; oscillatory associative memory is a funded hardware field. A register has no drift, holds state at zero power, and already has exponential range (binary is a modular code). The one open question is whether physical similarity search beats digital at large N, which needs hardware numbers, not simulation.

## Latest gate: Gate 4 — ping the memory with a goal, listen for the vector to it

![gate 4](results/gate4_summary.png)

Full write-up: **[GATE4.md](GATE4.md)**. Prompted by Sol's correction to Gate 3: an uncoupled oscillator's response operator *rotates* with its stored phase, so a fixed probe gets a phase-dependent answer with no coupling at all. Singular values (Gate 3's measure) are blind to that.

The task: path-integrate a position (Gate 2). Then ping the bank with the phase pattern a goal location *would* produce, and listen to each unit's phase advance. To first order that advance is sin(goal phase − stored phase), which interferes the goal against the memory.

- **It works without coupling (K1 passes).** Goal-vector error is 0.027 (box width 1), against 0.235 for answering "the goal is here".
- **It is nearly ideal (K2 narrowly fails).** Reading the phases directly gives 0.018, so the physical ping loses about a third.
- **Listening barely disturbs the memory.** A small ping adds 0.0003 to position error.
- **Vortex coupling adds nothing (K3 fails).** Coupling switched on by the ping ties at best; steady coupling is 2.4× worse. This is the fourth gate in which vortex coupling does not earn its place.
- **Correction to Gates 2 and 3:** with a capable (kNN) reader, the uncoupled bank's absolute position error is **0.021, not 0.198**. The linear reader understated the memory by about 9×. Each unit drifts 0.71 rad, yet 40 units together pin position: redundancy does the error correction for free.

**What survives:** integrate silently in uncoupled isochronous oscillators; query by stimulating and listening; no coupling needed. The ping's value is the interface. A reader that can only inject into a port and listen, never inspect internal state, can still ask "where is the goal from here?". That is the situation of a physical substrate. None of this is new neuroscience (oscillatory interference; goal vectors from grid codes, Bush et al. 2015).

## Gate 3 — memory vs control, and event-gated coupling

![gate 3](results/gate3_summary.png)

Full write-up: **[GATE3.md](GATE3.md)**, which now opens with Sol's correction. Coupling makes stored position change the *strength* of the bank's responses (rotation-invariant control) at the direct price of integration accuracy, and coupling switched on only 10% of the time gets 7–9× more of that control at the same error. **Corrected:** the original claim that an uncoupled bank's operator "cannot see what is stored" is wrong. It rotates with the stored phase, which Gate 4 shows is enough for a useful query. Gate 3's error axis also used a linear reader (see Gate 4).

- Rotation-invariant control is exactly 0 without coupling and grows linearly with κ; no coupling lowers integration error.
- A global rotation leaves singular values untouched under linear coupling (4×10⁻¹⁶) and changes them by 3.8% under vortex coupling.
- Not pre-registered: vortex beats linear on both axes at every κ ≥ 0.01.

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
pip install numpy scipy matplotlib
python vmn_checks.py         # ~70 s
python make_figure.py
python gate1_layer.py        # ~13 min
python gate1b_followups.py   # ~8 min
python make_gate1_figure.py
python gate2_path.py         # ~5 min
python make_gate2_figure.py
python gate3_tradeoff.py     # ~4 min
python make_gate3_figure.py
python gate4_ping.py         # ~1 min
python make_gate4_figure.py
```

Files: `vmn_core.py` (vortex dynamics, exact Jacobians, oscillators) · `vmn_checks.py` (maths checks) · `gate1_layer.py`, `gate1b_followups.py` (reservoir test) · `gate2_path.py` (path integration) · `gate3_tradeoff.py` (memory vs control) · `gate4_ping.py` (ping-and-listen) · `make_figure.py`, `make_gate1_figure.py` … `make_gate4_figure.py` · `results/`.
