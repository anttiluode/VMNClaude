"""Draw results/gate1_summary.png from the Gate 1 and Gate 1b receipts."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

g1 = json.load(open("results/gate1_receipt.json"))
g1b = json.load(open("results/gate1b_receipt.json"))
col = {"vortex": "#2a6fdb", "linear": "#d1495b", "none": "#888888", "esn": "#222222"}
fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))

# (a) memory capacity vs distance to measured onset, co-rotating units, beta = 0
for mode in ("vortex", "linear", "none"):
    rows = g1["sweep"][f"{mode}_beta0"]
    ax[0].plot([r["delta"] for r in rows], [r["MC"] for r in rows], "o-", color=col[mode], label=mode)
ax[0].axhline(g1["esn"]["test"]["MC"], color=col["esn"], ls="--", lw=1, label="ESN")
ax[0].axvline(0, color="#bbbbbb", lw=0.8)
ax[0].set_xlabel("delta = mu - measured onset")
ax[0].set_ylabel("memory capacity (test)")
ax[0].set_title("Gate 1: units all turning the same way")
ax[0].legend(frameon=False, fontsize=8)

# (b) co- vs counter-rotating, best over delta, beta = 0
perf = g1b["B_counter_rotation"]["performance"]
labels, mc, na, colors = [], [], [], []
for mode in ("vortex", "linear"):
    labels += [f"{mode}\nsame sense", f"{mode}\nmixed sense"]
    mc += [g1["summary"][f"{mode}_beta0"]["MC_best"], perf[f"{mode}_beta0_counter"]["MC_best"]]
    na += [g1["summary"][f"{mode}_beta0"]["NARMA_best"], perf[f"{mode}_beta0_counter"]["NARMA_best"]]
    colors += [col[mode], col[mode]]
x = np.arange(len(labels))
bars = ax[1].bar(x, mc, color=colors)
for i in (1, 3):
    bars[i].set_hatch("//")
    bars[i].set_edgecolor("white")
esn_best_mc = max(g1["esn"]["test"]["MC"], g1b["A_grid_edge"]["esn_rho0.9"]["MC_s0.03"])
ax[1].axhline(esn_best_mc, color=col["esn"], ls="--", lw=1, label="ESN memory capacity")
ax[1].set_xticks(x, labels, fontsize=8)
ax[1].set_ylabel("memory capacity (bars; test, best delta)")
ax[1].set_ylim(0, 48)
ax[1].set_title("Mixed rotation senses: vortex coupling wakes up")

ax2 = ax[1].twinx()
ax2.plot(x, na, "kD", ms=6, label="NARMA10 error")
ax2.axhline(g1["esn"]["test"]["NARMA"], color="#555555", ls=":", lw=1, label="ESN NARMA10 error")
ax2.set_ylabel("NARMA10 NRMSE (diamonds; lower is better)")
ax2.set_ylim(0.3, 0.6)
ax2.spines["top"].set_visible(False)
h1, l1 = ax[1].get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax[1].legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")

# (c) consistency: does the reservoir forget where it started?
for mode in ("vortex", "linear"):
    rows = g1["sweep"][f"{mode}_beta0"]
    c = np.array([max(r["CONS_MC_cfg"], 1e-17) for r in rows])
    ax[2].semilogy([r["delta"] for r in rows], c, "o-", color=col[mode], label=mode)
ax[2].axvline(0, color="#bbbbbb", lw=0.8)
ax[2].set_xlabel("delta = mu - measured onset")
ax[2].set_ylabel("state difference from two starts")
ax[2].set_title("Consistency breaks just above onset")
ax[2].legend(frameon=False, fontsize=8)

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("results/gate1_summary.png", dpi=130)
print("wrote results/gate1_summary.png")
