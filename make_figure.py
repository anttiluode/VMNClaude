"""Draw results/vmn_summary.png from results/receipt.json."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

d = json.load(open("results/receipt.json"))
ink, a1, a2 = "#222222", "#2a6fdb", "#d1495b"
fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))

# (a) event update rank vs N
t = d["C2_event_tidal_update"]["rank_vs_N_fixed_density"]
N = [r["N"] for r in t]
ax[0, 0].semilogx(N, [2 * n for n in N], ":", color="#999999", label="dimension 2N")
ax[0, 0].semilogx(N, [r["median_95pct_energy_rank"] for r in t], "o-", color=a1, label="update, 95% energy rank")
ax[0, 0].semilogx(N, [r["median_participation_rank"] for r in t], "s-", color=a2, label="update, participation rank")
ax[0, 0].set_yscale("log")
ax[0, 0].set_xlabel("number of vortices N (fixed density)")
ax[0, 0].set_title("C2  one new vortex: tidal update rank saturates")
ax[0, 0].legend(frameon=False, fontsize=8)

# (b) finite horizon
rows = [r for r in d["C3b_update_rank_vs_horizon"]["rows"] if r["kick"] == 0.25]
T = [r["T"] for r in rows]
ax[0, 1].plot(T, [r["propagator_participation_rank"] for r in rows], "o-", color="#999999", label="propagator itself")
ax[0, 1].plot(T, [r["update_participation_rank"] for r in rows], "s-", color=a2, label="change caused by the event")
ax[0, 1].set_yscale("log")
ax[0, 1].set_xlabel("horizon T")
ax[0, 1].set_ylabel("participation rank (of 128)")
ax[0, 1].set_title("C3  finite horizon, N = 64")
ax[0, 1].legend(frameon=False, fontsize=8)

# (c) vortex pair in strain
fp = d["C4_vortex_pair_as_neuron"]["free_pair"]["period_table"]
x = [np.log(1 / r["H_minus_Hs"]) for r in fp]
ax[1, 0].plot(x, [r["period"] for r in fp], "o-", color=a1, label="free pair (measured)")
s = d["C4_vortex_pair_as_neuron"]["free_pair"]["predicted_2/lambda"]
ax[1, 0].plot(x, [fp[-1]["period"] + s * (xi - x[-1]) for xi in x], "--", color=ink, lw=1, label="slope 2/lambda = 1/e")
ax[1, 0].set_xlabel("log 1/(H - H_saddle)")
ax[1, 0].set_ylabel("rotation period")
ax[1, 0].set_title("C4  co-rotating pair in strain: log onset")
ax[1, 0].legend(frameon=False, fontsize=8)

# (d) memory near Hopf
ab = d["C5_memory_near_hopf"]["above_onset"]
mu = [r["mu"] for r in ab]
ax[1, 1].loglog(mu, [abs(r["phase_shift"]) for r in ab], "o", color=a1, label="phase kept after kick (measured)")
ax[1, 1].loglog(mu, [abs(r["predicted_-beta*rho0/sqrt(mu)"]) for r in ab], "-", color=ink, lw=1, label="beta rho0 / sqrt(mu)")
ax[1, 1].set_xlabel("mu (distance above Hopf onset)")
ax[1, 1].set_title("C5  Stuart-Landau: memory diverges at onset")
ax[1, 1].legend(frameon=False, fontsize=8)

for a in ax.ravel():
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("results/vmn_summary.png", dpi=130)
print("wrote results/vmn_summary.png")
