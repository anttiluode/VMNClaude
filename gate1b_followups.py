"""Gate 1b - two follow-ups to Gate 1, run after seeing its result.

A. Grid edge: Gate 1 tuning chose the smallest input scale s = 0.1 almost everywhere.
   Re-test each variant (and the ESN) at s = 0.03 and 0.01, at its own best delta
   and kappa, so the verdicts are not an artefact of the grid edge.

B. Rotating-wave mechanism: antilinear coupling links z_j (turning at +w_j) to
   z_k (turning at +w_k) through conj(z_j), which turns at -w_j. With all units
   turning the same way that link is off-resonant by w_k + w_j, so it averages
   out; linear coupling is resonant when w_k ~ w_j. Prediction: mixing rotation
   senses (w_k -> sign(g_k) w_k, as vortices of both signs have) makes antilinear
   coupling resonant for opposite-sense pairs. Measured by the onset shift and by
   MC / NARMA, tuned on tune seeds and tested on sweep seeds exactly as in Gate 1.
"""
import json
import time

import numpy as np

import gate1_layer as g

t0 = time.time()
rec = json.load(open("results/gate1_receipt.json"))
TUNE, TEST = (0, 1), (10, 11, 12)
out = {"A_grid_edge": {}, "B_counter_rotation": {}}

# A ------------------------------------------------------------------------
for key, row in rec["summary"].items():
    if key == "esn":
        continue
    mode, b = key.split("_beta")
    beta = float(b)
    res = {}
    for task in ("MC", "NARMA"):
        k, _ = rec["tuning"][key]["chosen"][task]
        d = row[f"{task}_best_delta"]
        for s in (0.03, 0.01):
            vals = [g.eval_vmn(sd, mode, beta, k, s, d)[task] for sd in TEST]
            res[f"{task}_s{s}"] = float(np.mean(vals))
        res[f"{task}_s0.1_gate1"] = row[f"{task}_best"]
    out["A_grid_edge"][key] = res
    print("A", key, res, f"{time.time()-t0:.0f}s", flush=True)
esn = {}
for task in ("MC", "NARMA"):
    for s in (0.03, 0.01):
        esn[f"{task}_s{s}"] = float(np.mean([g.eval_esn(sd, 0.9, s)[task] for sd in TEST]))
    esn[f"{task}_s0.1_gate1"] = rec["esn"]["test"][task]
out["A_grid_edge"]["esn_rho0.9"] = esn
print("A esn", esn, flush=True)

# B ------------------------------------------------------------------------
onset = []
for kappa in (0.3, 1.0, 3.0):
    for counter in (False, True):
        for mode in ("vortex", "linear"):
            shifts = []
            for sd in TEST:
                rng = np.random.default_rng(1000 + sd)
                M, gg = g.vortex_matrix(rng)
                w = rng.uniform(0.5, 1.5, g.N)
                if counter:
                    w = w * np.sign(gg)
                shifts.append(g.onset_mu(M, w, kappa, mode))
            onset.append({"kappa": kappa, "counter_rotating": counter, "mode": mode,
                          "onset_shift_mu_c": float(np.mean(shifts))})
out["B_counter_rotation"]["onset"] = onset
for r in onset:
    print("B onset", r, flush=True)

perf = {}
for mode in ("vortex", "linear"):
    for beta in (0.0, 2.0):
        key = f"{mode}_beta{beta:g}_counter"
        grid = []
        for k in (0.3, 1.0, 3.0):
            for s in (0.03, 0.1):
                rs = [g.eval_vmn(sd, mode, beta, k, s, -0.05, counter=True) for sd in TUNE]
                grid.append({"kappa": k, "s": s,
                             "MC_val": float(np.mean([r["MC_val"] for r in rs])),
                             "NARMA_val": float(np.mean([r["NARMA_val"] for r in rs]))})
        cmc = max(grid, key=lambda r: r["MC_val"])
        cna = min(grid, key=lambda r: r["NARMA_val"])
        rows = []
        for d in (-0.2, -0.1, -0.05, -0.02, 0.0, 0.03):
            rmc = [g.eval_vmn(sd, mode, beta, cmc["kappa"], cmc["s"], d, counter=True) for sd in TEST]
            rna = [g.eval_vmn(sd, mode, beta, cna["kappa"], cna["s"], d, counter=True) for sd in TEST]
            rows.append({"delta": d, "MC": float(np.mean([r["MC"] for r in rmc])),
                         "NARMA": float(np.mean([r["NARMA"] for r in rna])),
                         "CONS_MC_cfg": float(np.nanmean([r["CONS"] for r in rmc])),
                         "mu_c": float(np.mean([r["mu_c"] for r in rmc]))})
        perf[key] = {"chosen_MC": [cmc["kappa"], cmc["s"]], "chosen_NARMA": [cna["kappa"], cna["s"]],
                     "rows": rows,
                     "MC_best": max(r["MC"] for r in rows), "NARMA_best": min(r["NARMA"] for r in rows)}
        print("B", key, perf[key]["MC_best"], perf[key]["NARMA_best"], f"{time.time()-t0:.0f}s", flush=True)
out["B_counter_rotation"]["performance"] = perf
out["seconds"] = round(time.time() - t0, 1)
json.dump(out, open("results/gate1b_receipt.json", "w"), indent=2)
print("done", out["seconds"])
