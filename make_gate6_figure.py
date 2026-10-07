import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open("results/gate6_receipt.json"))
S, AMPS = R["summary_mean_sd"], R["setup"]["amps"]
B = json.load(open("results/gate6b_receipt.json"))
col = {"single": "#8a3a3a", "repeat": "#d9822b", "sign_flip": "#2a6fdb", "conjugate": "#3a8a5c"}
lab = {"single": "single ping", "repeat": "ping twice, same sign (no undo)",
       "sign_flip": "ping, listen, ping with -sign (undo)", "conjugate": "ping, listen, conjugate kick (undo)"}

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
for s in col:
    e = [S[f"uniform|{s}|{a:g}"]["answer_error"][0] for a in AMPS]
    d = [S[f"uniform|{s}|{a:g}"]["phase_disturbance_rms"][0] for a in AMPS]
    ax[0].plot(d, e, "o-", color=col[s], label=lab[s])
    for a, x, y in zip(AMPS, d, e):
        if s == "single":
            ax[0].annotate(f"a={a:g}", (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7, color=col[s])
ax[0].axhline(S["oracle_uniform"], color="k", ls=":", lw=1, label=f"read phases directly ({S['oracle_uniform']:.3f})")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("memory damage left behind: RMS phase disturbance per unit (rad)")
ax[0].set_ylabel("goal-vector answer error (box = 1)")
ax[0].set_title("Same answer, much less damage: lower-left is better", fontsize=10)
ax[0].legend(fontsize=8, loc="upper right")

keys = [("real memory, noise on", "sigma0.05_window16_a0.3"), ("real memory, noise off", "sigma0_window64_a0.3"),
        ("on limit cycle, noise off", "on_limit_cycle_sigma0_window64_a0.3")]
sch = [("single", "#8a3a3a"), ("sign_flip", "#2a6fdb"), ("conjugate_current", "#7fbf95"), ("conjugate_prequery", "#3a8a5c")]
w = 0.2
for i, (lbl, k) in enumerate(keys):
    for j, (s, c) in enumerate(sch):
        ax[1].bar(i + (j - 1.5) * w, B[k][s], w, color=c, label=s.replace("_", " ") if i == 0 else None)
ax[1].set_yscale("log"); ax[1].set_xticks(range(3)); ax[1].set_xticklabels([k[0] for k in keys], fontsize=8)
ax[1].set_ylabel("residual phase disturbance (rad), a = 0.3")
ax[1].set_title("Exact undo needs the pre-query phase and a clean amplitude", fontsize=10)
ax[1].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4)
fig.suptitle("Gate 6: read the goal, listen, then restore the memory", fontsize=11)
fig.tight_layout(); fig.savefig("results/gate6_summary.png", dpi=130)
