import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = json.load(open("results/gate4_receipt.json"))["summary_mean_sd"]
rows = [("uncoupled ping, a=0.05", "none_a0.05"), ("uncoupled ping, a=0.3", "none_a0.3"),
        ("ping-gated vortex 0.1", "gated_vortex_0.1_a0.3"), ("ping-gated vortex 0.3", "gated_vortex_0.3_a0.3"),
        ("ping-gated linear 0.1", "gated_linear_0.1_a0.3"), ("ping-gated linear 0.3", "gated_linear_0.3_a0.3"),
        ("steady vortex 0.03", "fixed_vortex_a0.3")]
err = [S[k]["best"][0] for _, k in rows]
sd = [S[k]["best"][1] for _, k in rows]
dmg = [S[k]["position_error_after_knn"][0] - S[k]["position_error_before_knn"][0] for _, k in rows]
col = ["#3a8a5c", "#3a8a5c", "#2a6fdb", "#2a6fdb", "#d9822b", "#d9822b", "#8a3a3a"]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
y = range(len(rows))[::-1]
ax[0].barh(list(y), err, xerr=sd, color=col)
ax[0].axvline(S["predict_zero_error"], color="grey", ls="--", lw=1, label=f"answer 'goal is here' ({S['predict_zero_error']:.3f})")
orc = min(v[0] for v in S["oracle_int_none_0"].values())
ax[0].axvline(orc, color="k", ls=":", lw=1.2, label=f"oracle: read phases directly ({orc:.3f})")
ax[0].set_yticks(list(y)); ax[0].set_yticklabels([r[0] for r in rows])
ax[0].set_xscale("log"); ax[0].set_xlabel("goal-vector error (box width = 1), best reader, log scale")
ax[0].set_title("Ping the memory with a goal, listen for the vector to it"); ax[0].legend(fontsize=8, loc="lower right")
ax[1].barh(list(y), dmg, color=col)
ax[1].set_yticks(list(y)); ax[1].set_yticklabels([])
ax[1].axvline(0, color="k", lw=0.8)
ax[1].set_xlabel("read damage: position error added by the ping")
ax[1].set_title("Does listening disturb the memory?")
fig.tight_layout(); fig.savefig("results/gate4_summary.png", dpi=130)
