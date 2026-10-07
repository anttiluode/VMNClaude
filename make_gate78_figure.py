import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G7 = json.load(open("results/gate7_receipt.json"))["summary_mean_sd"]["a0.2"]
R8 = json.load(open("results/gate8_receipt.json"))
G8, zero = R8["summary_mean_sd"]["a0.2"], R8["summary_mean_sd"]["predict_zero_error"]
G8b = json.load(open("results/gate8b_receipt.json"))["summary_mean_sd"]

fig, ax = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1, 1.35]})
rows7 = [("no undo", "none", "#8a3a3a"), ("undo from the 2-number answer", "answer", "#d9822b"),
         ("undo from the TRUE answer", "answer_true", "#e8b27a"), ("sign-flip undo (memoryless)", "sign_flip", "#2a6fdb"),
         ("undo from per-unit pre-query copy", "prequery", "#3a8a5c")]
y = list(range(len(rows7)))[::-1]
ax[0].barh(y, [G7[k][0] for _, k, _ in rows7], xerr=[G7[k][1] for _, k, _ in rows7], color=[c for *_, c in rows7])
ax[0].set_yticks(y); ax[0].set_yticklabels([r[0] for r in rows7], fontsize=9)
ax[0].set_xlabel("lasting phase damage (rad), a = 0.2")
ax[0].set_title("Gate 7: the answer can't tell you how to undo the read", fontsize=10)

rows8 = [("full listener: every unit (4 t.u.)", G8["full"][0], "#3a8a5c"),
         ("1 summed channel", G8["chan_1"][0], "#8a3a3a"), ("2 summed channels", G8["chan_2"][0], "#8a3a3a"),
         ("4 summed channels", G8["chan_4"][0], "#8a3a3a"),
         ("one FDM wire, no noise (12 t.u. windows)", G8["wire_fdm_sigma0"][0], "#2a6fdb"),
         ("one FDM wire, noise 0.5", G8["wire_fdm_sigma0.5"][0], "#2a6fdb"),
         ("one FDM wire, noise 2", G8["wire_fdm_sigma2"][0], "#2a6fdb"),
         ("perfect demodulator, 12 t.u. windows (8b)", G8b["avg_long"][0], "#7fa8e8"),
         ("perfect demodulator, short windows (8b)", G8b["avg_short"][0], "#7fbf95")]
y = list(range(len(rows8)))[::-1]
ax[1].barh(y, [r[1] for r in rows8], color=[r[2] for r in rows8])
ax[1].axvline(zero, color="grey", ls="--", lw=1, label=f"answer 'goal is here' ({zero:.3f})")
ax[1].set_yticks(y); ax[1].set_yticklabels([r[0] for r in rows8], fontsize=9)
ax[1].set_xscale("log"); ax[1].set_xlabel("goal-vector answer error, a = 0.2 (log)")
ax[1].set_title("Gate 8: a cheap interface needs frequency multiplexing AND bandwidth", fontsize=10)
ax[1].legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig("results/gate78_summary.png", dpi=130)
