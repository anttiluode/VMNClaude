# Vortex – Matrix – Neuron

**Working note, 7 October 2026. Claude (Opus 5.5).**

The question: is there exact mathematics linking a vortex, a matrix and a neuron, with the nonlinearity doing the real work? This note gives four links. Each one is checked numerically by `vmn_checks.py`, and every number below comes from `results/receipt.json`.

Ledger first:

| | Status |
|---|---|
| Point-vortex equations, Hamiltonian structure, Adler equation, SNIC vs homoclinic onset, Stuart–Landau, Landau equation for the cylinder wake | **Classical.** Cited, not claimed. |
| Antilinear complex-symmetric form of the vortex response (§1) | Derived here. Standard algebra; probably known in pieces. Verified to machine precision. |
| Tidal update law and N-independent rank of an event (§2) | Derived here. Exact one-step result, verified. Finite-horizon behaviour is **measured**, not proved. |
| Co-rotating pair in strain as a neuron, with two onset classes (§3) | Derivation is elementary; verified numerically. |
| Memory scaling β/√μ near the Hopf point (§4) | Standard phase-reduction result, verified. Its use as an **explanation** of the Kármán computing peak is a hypothesis, untested on that system. |

![summary](results/vmn_summary.png)

---

## 0. Why the nonlinearity has to be there

For any dynamics $`\dot x=f(x)`$ the response to the next small input is the Jacobian $`K(x)=Df(x)`$. An event that moves the state by $`a`$ changes it by

```math
\Delta K \;=\; Df(x+a)-Df(x)\;=\;D^2f(x)[a]+O(|a|^2).
```

If $`f`$ is linear, $`D^2f=0`$ and **no event can ever change how the system responds**. A linear system can carry a signal, but it cannot store an experience in its response geometry. The second derivative of the vector field is exactly the channel through which the past rewrites the operator. Kompressori's ΔJ, a neuron's phase response curve, and a vortex reshaping the flow are all this one term.

---

## 1. Vortex → Matrix: the response is an antilinear complex-symmetric map

Point vortices $`z_k=x_k+iy_k`$ with circulations $`\Gamma_k`$ obey (Helmholtz–Kirchhoff)

```math
\dot{\bar z}_k=\frac{1}{2\pi i}\sum_{j\ne k}\frac{\Gamma_j}{z_k-z_j}.
```

Differentiating gives the tangent dynamics in complex form:

```math
\delta\dot z \;=\; M\,\overline{\delta z},\qquad
M_{kj}=\frac{i}{2\pi}\frac{\Gamma_j}{(\bar z_k-\bar z_j)^2}\;(k\ne j),\qquad
M_{kk}=-\sum_{j\ne k}M_{kj}.
```

Three facts follow, and all three are checked:

1. **The response is antilinear.** A perturbation $`\delta z`$ produces a velocity proportional to its *conjugate*. Geometrically, each 2×2 block is a pure strain (symmetric, traceless), never a rotation. This is because $`\log|z|`$ is harmonic. Measured: $`\operatorname{tr}K=0`$, and the antilinear form matches the real Jacobian to $`3\times10^{-16}`$.
2. **Weighted by circulation, it is complex symmetric:** $`\operatorname{diag}(\Gamma)M=(\operatorname{diag}(\Gamma)M)^{\mathsf T}`$, measured to $`10^{-16}`$. Antilinear maps of this type are diagonalised by the Takagi factorisation, not ordinary eigenvectors.
3. **Frequencies are square roots of a linear spectrum.** Applying the map twice gives a complex-linear map, so

```math
\operatorname{spec}(K)^2=\operatorname{spec}(M\bar M)\cup\overline{\operatorname{spec}(M\bar M)}.
```

Verified exactly (to 9 digits) on random configurations. The oscillation frequencies of any vortex arrangement come from one complex matrix, $`M\bar M`$.

**What this buys:** a vortex system *is* a matrix system, specifically a strain-only, antilinear one. Any "vortex–matrix coupling" design should keep that structure; an ordinary symmetric or rotation coupling would be a different physics.

---

## 2. An event writes a tidal update, and its rank does not grow with N

Add one vortex of strength $`g`$ at position $`w`$ (the "event"). The response among the existing vortices changes only in the diagonal blocks, because the new vortex enters each old vortex's velocity separately:

```math
\Delta K=\operatorname{blockdiag}\big(B_1,\dots,B_N\big),\qquad
\sigma_1(B_k)=\sigma_2(B_k)=\frac{|g|}{2\pi\,|z_k-w|^2}.
```

Each block is a tidal strain with both singular values $`|g|/2\pi r_k^2`$. Two consequences:

- **The update a vortex feels does not depend on its own circulation.** Measured: off-block-diagonal norm exactly 0; singular values match the formula to $`10^{-14}`$.
- **The effective rank saturates.** The participation rank is $`2\big(\sum r_k^{-4}\big)^2/\sum r_k^{-8}`$. In two dimensions $`\sum r^{-4}`$ converges at fixed density, so the rank is set by the nearest neighbours, not by N.

Measured at fixed density, median over 40 draws:

| N | dimension | participation rank | 95%-energy rank |
|---:|---:|---:|---:|
| 25 | 50 | 6.1 | 16 |
| 400 | 800 | 5.8 | 24.5 |
| 6400 | 12800 | 5.7 | 25.5 |

The 1/r² tidal law is why a local event is cheap to describe. This is the vortex version of the one-step bound in the Kompressori note (one cell touches at most five Jacobian rows), but here the cause is a decay law rather than a finite stencil.

### Finite horizon (measured, not proved)

The one-step result is exact. Over a horizon T the effect spreads, so I integrated the variational equation for N smoothed vortices (δ = 0.05, density 1, minimum spacing 0.45), kicked the central vortex, and compared propagators. Medians over 4 configurations:

| N | dimension | propagator 95% rank | update 95% rank | update participation rank |
|---:|---:|---:|---:|---:|
| 16 | 32 | 25 | 5 | 2.5 |
| 36 | 72 | 55.5 | 7 | 3.7 |
| 64 | 128 | 93 | 4 | 2.3 |
| 100 | 200 | 144.5 | 4.5 | 2.3 |

At T = 1 the propagator is nearly full rank and the update stays at 4–7 directions, independent of N. These magnitudes are close to Kompressori's field (95% rank 5 of 100, participation 2.75), in a completely different system.

**Important caveat, from the horizon sweep at N = 64.** As T grows, the propagator itself collapses: its participation rank falls 96 → 51 → 23 → 5.6 for T = 0.5, 1, 2, 4, because chaotic stretching aligns everything with a few unstable directions. The update stays at 1.5–2.9. So the comparison "low-rank change of a high-rank operator" is meaningful only for horizons shorter than the stretching time. At long horizons, *everything* is low rank and the result says nothing special. The same caveat should be checked in Kompressori.

Tiny kicks (0.01) and finite kicks (0.25) give nearly the same rank, so the low rank is not a large-amplitude effect.

---

## 3. Vortex → Neuron: a co-rotating pair in strain is an excitable oscillator

Two vortices with total circulation $`\Gamma`$, placed in a pure strain $`e`$ (an external flow $`u=ex,\ v=-ey`$). Their separation $`z=re^{i\varphi}`$ obeys exactly

```math
\dot z=\frac{i\Gamma}{2\pi\bar z}+e\,\bar z
\quad\Longleftrightarrow\quad
\dot\varphi=\frac{\Gamma}{2\pi r^2}-e\sin2\varphi,\qquad \dot r=e\,r\cos2\varphi .
```

**Frozen separation: the Adler equation.** At fixed r this is the Adler / theta-neuron phase equation. Its "drive" is the self-rotation $`a=\Gamma/2\pi r^2`$, and its "inhibition" is the strain. Below threshold ($`a<e`$) the pair locks: the neuron rests. Above, it rotates: the neuron fires, with period $`\pi/\sqrt{a^2-e^2}`$. That is the square-root onset of a type-I (SNIC) neuron. Measured periods match the formula to 4–5 digits, down to 0.2% above threshold.

**Free separation: a different neuron class.** Letting r move, the system is Hamiltonian with $`H=-\frac{\Gamma}{4\pi}\log r^2+e\,xy`$. Two saddles sit at $`r_s=\sqrt{\Gamma/2\pi e}`$ with eigenvalues $`\pm2e`$. Bound orbits encircle the origin; past the separatrix the pair is torn apart. Near the separatrix the period grows **logarithmically**:

```math
T\approx\frac{2}{\lambda}\log\frac{1}{H-H_s}+\text{const},\qquad \lambda=2e .
```

Measured slope 20.03 against predicted 20.0 (e = 0.05). On the far side of the separatrix, the pair separates without bound.

So one vortex pair contains both standard routes to firing. The frozen version is SNIC (square-root); the full version is a heteroclinic loop (logarithmic). Which one a "vortex neuron" shows depends on whether its separation is pinned, and that is a design choice.

---

## 4. Neuron → Matrix: memory of a kick, and why the Hopf point is special

The amplitude of the Kármán vortex street just past onset obeys the Landau equation (Provansal, Mathis & Boyer, *J. Fluid Mech.* 182, 1987). That is the Stuart–Landau oscillator

```math
\dot z=(\mu+i\omega)z-(1+i\beta)|z|^2z,
```

where μ is the distance past the Hopf point and β is the **shear**: how much the frequency depends on amplitude. So the vortex street, reduced to its amplitude, is literally an oscillator neuron.

Kick the amplitude by $`\rho_0`$. The amplitude relaxes in time $`1/2\mu`$, but because of shear the phase does not come back:

```math
\Delta\varphi_\infty=-\frac{\beta\,\rho_0}{\sqrt{\mu}} .
```

Measured with the same kick at μ = 0.4, 0.1, 0.025, 0.00625: the agreement is better than 0.01%, and the stored phase grows 8× as μ shrinks 64×. Below onset (μ < 0) a kick decays at rate |μ|, also measured exactly. **Both the memory time and the stored phase diverge at the Hopf point, from both sides.**

This is a candidate explanation, not a test, for the result in *Computing with vortices* (arXiv:2001.08502): memory and nonlinear processing of a cylinder wake peak near the onset of shedding. That onset is this Hopf point. Testing it would mean checking whether the wake's measured memory curve follows 1/|μ|.

**The Hamiltonian limit is a different kind of memory.** An ideal vortex has no damping: μ → 0 with nothing pulling the amplitude back. A satellite vortex orbiting circulation Γ has $`\omega=\Gamma/2\pi r^2`$, so a radial kick δr changes its frequency permanently, and the phase error grows **linearly forever**:

```math
\Delta\varphi(t)=-\frac{2\Gamma\,\delta r}{2\pi r^3}\,t .
```

Measured slope matches to 0.015%. A dissipative neuron stores a kick as a finite phase shift. A Hamiltonian vortex stores it as a frequency, which is an integrator. "Memory that doesn't fade" is exact here, but so is the price: errors accumulate too.

---

## 5. What the three have in common

The nonlinearity Antti sensed is the same term in all three, and it has a name: **shear**, frequency depending on amplitude.

| | Shear term |
|---|---|
| Vortex | self-rotation $`\Gamma/2\pi r^2`$ depends on separation |
| Neuron (Stuart–Landau) | β, phase speed $`\omega-\beta r^2`$ |
| Kompressori field | amplitude-dependent wave speed $`c_0^2/(1+\alpha\phi^2)`$ and the cubic potential (Duffing) |

The Kompressori row is an observation about its equations; whether switching those terms off removes its susceptibility windows is untested.

Shear is what turns an amplitude event into a lasting change of timing (§4) and what makes the response operator state-dependent (§0). Tidal 1/r² decay is what keeps the change local and low-rank (§2). Antilinear strain coupling is what vortex interaction actually looks like as a matrix (§1).

**On "the vortex is an extra step":** the checks support that intuition. Everything computationally relevant here reduces to complex oscillators with shear, coupled antilinearly through a 1/r² strain kernel. A fluid is one way to get that. It is not the only one, and simulating a fluid to get it is overhead.

An aside, flagged as an identification rather than a result: in a fluid, the local rotation rate of an element is exactly half the vorticity. The vorticity field is literally a local clock-rate field, the object the Clockfield line kept reaching for.

---

## 6. Next gate (proposed, not done)

A VMN unit layer:

```math
\dot z_k=(\mu+i\omega_k)z_k-(1+i\beta)|z_k|^2z_k+\sum_{j\ne k} \frac{i\,g_j}{2\pi\,(\bar p_k-\bar p_j)^{2}}\;\overline{z_j}+u_k(t),
```

with units at fixed 2D positions $`p_k`$, inputs $`u_k`$, and a learned linear readout. The coupling is the off-diagonal part of $`M`$ from §1, with the units' oscillator states standing in for vortex displacements. That substitution is the design step being tested, not a derivation.

Kill conditions, fixed in advance:

1. Against the same layer with ordinary linear coupling ($`z_j`$ instead of $`\bar z_j`$): if antilinear coupling gives no gain on memory-capacity or NARMA tasks, the vortex structure adds nothing.
2. Against LinOSS and coRNN at equal parameter count on long-memory tasks: no gain → not useful as an AI layer.
3. Sweep μ through zero: if task memory does not peak near μ = 0 as §4 predicts, §4's explanation is wrong for this system.

## Run

```bash
pip install numpy matplotlib
python vmn_checks.py      # ~70 s, writes results/receipt.json
python make_figure.py     # writes results/vmn_summary.png
```
