"""Draw results/gate2_summary.png from results/gate2_receipt.json."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

r = json.load(open("results/gate2_receipt.json"))
V = r["vmn"]
t = np.arange(0, r["setup"]["steps"] + 1, 20)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for key, lab, c, ls in (("w00_none_beta0", "oscillators, no coupling", "#888888", "-"),
                        ("w00_vortex_beta0", "vortex coupling", "#2a6fdb", "-"),
                        ("w00_linear_beta0", "linear coupling", "#d1495b", "-"),
                        ("w00_none_beta2", "no coupling, shear beta = 2", "#888888", ":")):
    ax[0].plot(t, V[key]["curve"], ls, color=c, marker="o", ms=3, label=lab)
ax[0].plot(t, r["esn"]["curve"], "--", color="#222222", label="echo state network")
ax[0].axhline(0.3826, color="#bbbbbb", lw=0.8)
ax[0].text(2, 0.388, "always guessing the centre", fontsize=7, color="#888888")
ax[0].set_xlabel("steps since the start landmark")
ax[0].set_ylabel("position error (box = 1)")
ax[0].set_title("Path integration, both rotation senses (w0 = 0)")
ax[0].legend(frameon=False, fontsize=8)

labels, vals, cols = [], [], []
for w0 in (0, 3):
    for m, c in (("none", "#888888"), ("vortex", "#2a6fdb"), ("linear", "#d1495b")):
        labels.append(f"{m}\nw0={w0}")
        vals.append(V[f"w0{w0}_{m}_beta0"]["test_end_error"])
        cols.append(c)
x = np.arange(len(labels))
ax[1].bar(x, vals, color=cols)
ax[1].set_xticks(x, labels, fontsize=8)
ax[1].set_ylim(0.15, 0.22)
ax[1].set_ylabel("end-of-episode error (test)")
ax[1].set_title("Coupling adds nothing beyond seed noise (beta = 0)")
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("results/gate2_summary.png", dpi=130)
print("wrote results/gate2_summary.png")
