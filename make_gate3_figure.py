import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open("results/gate3_receipt.json"))
col = {"vortex": "#2a6fdb", "linear": "#d9822b"}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for m in ("vortex", "linear"):
    f = R["fixed"][m]; g = R["gated"][m]
    ax[0].plot([q["control"] for q in f[1:]], [q["error"] for q in f[1:]], "-o", color=col[m], label=f"{m}, fixed κ")
    ax[0].plot([q["control"] for q in g], [q["error"] for q in g], "s", color=col[m], mfc="white", ms=8, mew=2,
               label=f"{m}, event-gated (10% duty)")
    ks = [q["kappa"] for q in f[1:]]
    ax[1].loglog(ks, [q["control"] for q in f[1:]], "-o", color=col[m], label=f"{m}: control")
ax[0].axhline(R["fixed"]["vortex"][0]["error"], color="k", ls=":", lw=1, label="uncoupled (control = 0)")
ax[0].axhline(0.3826, color="grey", ls="--", lw=1, label="guess box centre")
ax[0].set_xscale("log"); ax[0].set_xlabel("control: how much stored position changes the response operator")
ax[0].set_ylabel("path-integration end error"); ax[0].set_title("Memory vs control: down-and-right is better")
ax[0].legend(fontsize=7.5)
ax2 = ax[1].twinx()
for m in ("vortex", "linear"):
    f = R["fixed"][m]
    ax2.semilogx([q["kappa"] for q in f[1:]], [q["error"] for q in f[1:]], "--", color=col[m], alpha=0.7)
ax2.set_ylabel("end error (dashed)")
ax[1].set_xlabel("coupling κ"); ax[1].set_ylabel("control (solid)"); ax[1].set_title("Both rise with κ: the trade-off")
ax[1].legend(fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig("results/gate3_summary.png", dpi=130)
