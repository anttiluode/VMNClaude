import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = json.load(open("results/gate5_receipt.json"))["partA"]
B = json.load(open("results/gate5b_receipt.json"))
NS = ["1", "2", "4", "8", "16", "32"]
col = {"same": "#8a3a3a", "left_right": "#d9822b", "random": "#2a6fdb", "antipodal": "#3a8a5c"}
lab = {"same": "same direction", "left_right": "left-right alternating (brain)",
       "random": "random direction", "antipodal": "antipodal (+D, -D)"}

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw={"width_ratios": [1.1, 1.1, 1.2]})
for i, (key, title) in enumerate([("beta0_a0.1_int8", "stationary, weak ping (a = 0.1), no shear"),
                                  ("beta0.5_a0.1_int8", "same, with shear (beta = 0.5)")]):
    for p in col:
        y = [A[key][p][n]["mean_abs_shift"][0] for n in NS]
        ax[i].plot([int(n) for n in NS], [max(v, 1e-5) for v in y], "o-", color=col[p], label=lab[p])
    ax[i].set_xscale("log", base=2); ax[i].set_yscale("log")
    ax[i].set_xlabel("number of sweeps (pings)"); ax[i].set_title(title, fontsize=10)
    ax[i].set_ylim(5e-5, 0.1)
ax[0].set_ylabel("read damage: shift of the stored position")
ax[0].legend(fontsize=8)

rows = [("same direction", "same"), ("left-right (brain)", "left_right"), ("random", "random"),
        ("antipodal", "antipodal"), ("left-right, pair shares heading", "left_right_paired"),
        ("antipodal, pair shares heading", "antipodal_paired")]
G5 = json.load(open("results/gate5_receipt.json"))["partB"]["a0.1"]
vals = [(B["a0.1"][k][0] if k in B["a0.1"] else G5[k]["mean_abs_shift"][0]) for _, k in rows]
c2 = ["#8a3a3a", "#d9822b", "#2a6fdb", "#7fbf95", "#e8b27a", "#3a8a5c"]
y = list(range(len(rows)))[::-1]
ax[2].barh(y, vals, color=c2)
ax[2].set_yticks(y); ax[2].set_yticklabels([r[0] for r in rows], fontsize=9)
ax[2].set_xscale("log"); ax[2].set_xlabel("read damage after 100 sweeps while navigating (a = 0.1)")
ax[2].set_title("navigating: heading turns ~37 deg between sweeps", fontsize=10)
fig.suptitle("Gate 5: does left-right alternation protect the phase memory from its own reads?  No.", fontsize=11)
fig.tight_layout(); fig.savefig("results/gate5_summary.png", dpi=130)
