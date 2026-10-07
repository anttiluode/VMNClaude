# Who needs ping-queryable physical memory? An honest market check

**7 October 2026. Claude (Opus 5.5).** Follows [Gate 4](GATE4.md). Question from Antti: is there anything that needs "integrate silently in oscillators, query by pinging" more than it needs an ordinary digital register?

**Short answer: almost nothing, and where it is needed, it already exists.** The principle is real and deployed, but in mature engineering niches where we would be late, small entrants. This note shows the reasoning so the verdict can be checked, and names the one test that could change it.

## The three things a register does better

1. **No drift.** A digital integrator adds no error of its own; every bit of error comes from the sensor. The oscillator bank in Gates 2–4 *adds* error, because its phases diffuse (0.71 rad RMS per unit by the end of an episode). Redundancy across 40 units cleans most of it up, but a register starts from zero added error and needs no redundancy at all.
2. **Holds state at zero power.** A limit cycle needs energy every moment it oscillates (μ > 0 is a pump). Flash, FRAM or MRAM hold a value with the power off. Any application that stores something for a long time and reads it rarely, which is where pinging would be most useful, is exactly where a sustained oscillator loses.
3. **Range is already exponential.** The grid code's famous property, exponential range from a few modules, is something binary already has: a binary counter *is* a modular code with moduli 2, 4, 8, … So "exponential range from bounded state" is not an advantage over digital.

So the oscillator has to win on something else: the input is *already* a frequency, the reader *cannot* be wired to the state, or the *comparison* (not the storage) is the expensive part. Each of these has an existing industry.

## Where the principle is already used

| Our ingredient | Where it already lives | What that means for us |
|---|---|---|
| Phase of a voltage-controlled oscillator as an integrator | **VCO-based ADCs.** The input sets the oscillator's frequency, the phase accumulates it, and counting phase gives first-order noise shaping for free ([Straayer, MIT thesis](https://dspace.mit.edu/handle/1721.1/47755); [Perrott's tutorial](https://bibliotheek.ehb.be:3124/r5/central_texas/cas_ssc/meetings/2012/101112/vco_based_quant_perrott.pdf); recent hybrid designs from the [ASU mixed-signal group](https://labs.engineering.asu.edu/mixedsignals/wp-content/uploads/sites/58/2023/10/NS-SAR-VCO-ESSCIRC_2023.pdf)) | Gate 2's "velocity sets frequency, phase stores the integral" is this circuit. It's mature, with a decade-plus of silicon. |
| Ping a physical element, listen to the echo for its state, no wires or battery | **Passive wireless SAW sensors and tags** interrogated by RF pulses, including shaft-angle and RPM sensing on rotating parts ([US7065459](https://patents.google.com/patent/US7065459), [SAW tachometer US9817014](https://patents.google.com/patent/US9817014)) | A rotating shaft *is* a phase integrator of angular velocity, and reading it by ping already exists. |
| One broadcast ping, many oscillators at different frequencies, read on a single line | **Frequency-multiplexed resonator readout,** used for example for kilopixel superconducting detector arrays in astronomy (from my own knowledge; not checked in this session) | Gate 4's per-unit phase advance, read on one wire, is this readout scheme. |
| Ping with a pattern, the substrate answers with similarity | **Oscillatory neural networks / oscillator associative memory,** an active hardware research field with dedicated reviews and digital/analog ONN chips ([Todri-Sanial et al.](https://research.tue.nl/en/publications/building-oscillatory-neural-networks-ai-applications-and-physical/); [fully connected digital ONN scaling](https://arxiv.org/abs/2504.20680); [Frontiers 2024](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2024.1307525/pdf)) | Gate 4's goal query is a tiny case of this, and the field is crowded with funded hardware teams. |
| Grid-cell path integration in hardware | **Neuromorphic implementations** exist, for example on spiking hardware for robots ([Kreiser et al.](https://blogs.ed.ac.uk/rai-nr/wp-content/uploads/sites/406/2019/07/Kreiser_RobustAI.pdf); [Frontiers 2020](https://frontiersin.org/article/10.3389/fnins.2020.00551/full)) | Novelty as a navigation system: none. |

## The decision table

| Situation | Winner | Why |
|---|---|---|
| Store a number, read it later | Register | No drift, zero hold power |
| Integrate a sensor that outputs a voltage | Register after an ADC, or a VCO-ADC | Both are commodity parts |
| Integrate a sensor that already outputs a frequency (resonant MEMS, quartz) | **The sensor's own phase, read by a counter** | The oscillator memory already exists for free; no bank needed |
| Read state through a wall or off a moving part, with no battery | **Passive SAW-style ping** | Already a product category |
| Read thousands of integrators on one wire | **Frequency-multiplexed resonators** | Already used in detector arrays |
| Compare a query against millions of stored patterns at once | *Open:* oscillator associative memory vs digital CAM/GPU | The only place the physics might beat digital, and an active research race |
| Understand how the brain navigates | The oscillator model | Scientific value, not a market |

## The one test that could change the verdict

The only unsettled row is the comparison step: "ping with a query, the substrate computes similarity to everything it holds, in one shot." Digital content-addressable memory is notoriously power-hungry, so a physical bank could win there. But it only wins with a hardware cost model: energy per query and area per stored item, against a digital baseline, as the number of stored items N grows. Python simulations can't measure that. The honest version needs component numbers from real oscillator hardware (published ONN chips report them) and a crossover N where physics beats a GPU nearest-neighbour search. If the crossover N is beyond any realistic workload, the line is closed for money.

## Recommendation

- **As a product, I'd stop here.** Every place ping-queryable oscillator memory is needed is already served, and the one open race (oscillator associative memory) has funded hardware teams, not a solo build edge.
- **What VMN did produce is worth keeping:** a clean measured architecture (integrate silently, query by ping, redundancy instead of coupling), four honest gates showing vortex coupling doesn't earn a place, and a reader correction that changed every absolute number by 9×.
- If you want one more step with a real chance of value, it is the crossover estimate above, done as a desk calculation from published chip numbers, before writing any code.
