"""plot_results.py - renders every figure from CSV files written by the Octave/MATLAB
pipeline (results/). No numbers are generated here; all data come from executed runs.
(The MATLAB equivalent is m19_plot_results.m; Octave plotting was unavailable headless.)"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, os
R = "results"; F = "figures_v2"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.dpi": 110, "savefig.dpi": 200, "axes.grid": True, "grid.alpha": .3})
rd = lambda f: pd.read_csv(os.path.join(R, f))
E = rd("st02_engineering_boit.csv")
def save(fig, n): fig.tight_layout(); fig.savefig(os.path.join(F, n), bbox_inches="tight"); plt.close(fig)

# F01 Brayton surfaces at eps_r = 0.95
S = rd("st09_surface_brayton.csv"); nr, nt = 111, 121
g = lambda c: S[c].values.reshape(nt, nr).T if False else S[c].values.reshape((nr, nt), order="F")
rp, T3 = g("rp"), g("T3"); feas = g("cv") <= 0
panels = [("P", r"$P=\eta_{th}$"), ("BT", r"$B_T=1-T_0/T_3$"), ("BA", r"$B_C$ (heat-exergy limit)"),
          ("AT", r"$A_T=\eta/B_T$"), ("AA", r"$A_C=\eta/B_C=\eta_{II}$"), ("GA", r"$G_C=1-A_C$")]
fig, ax = plt.subplots(2, 3, figsize=(11, 6.2))
xP = E.loc[0, ["xP_1", "xP_2"]].values; xAT = E.loc[0, ["xA_1", "xA_2"]].values; xAA = E.loc[1, ["xA_1", "xA_2"]].values
for a, (c, t) in zip(ax.ravel(), panels):
    Z = 1 - g("AA") if c == "GA" else g(c)
    cs = a.contourf(rp, T3, Z, 20, cmap="viridis"); fig.colorbar(cs, ax=a)
    a.contourf(rp, T3, (~feas).astype(float), levels=[.5, 1.5], colors="none", hatches=["///"])
    a.contour(rp, T3, g("T4"), levels=[1000], colors="w", linewidths=1)
    a.plot(*xP, "r*", ms=11, label=r"$x^*_P$"); a.plot(*xAT, "ws", ms=6, mec="k", label=r"$x^*_{A_T}$"); a.plot(*xAA, "c^", ms=6, mec="k", label=r"$x^*_{A_C}$")
    a.set_title(t); a.set_xlabel("$r_p$"); a.set_ylabel("$T_3$ (K)")
ax[0, 0].legend(fontsize=7, loc="lower right")
fig.suptitle(r"Brayton, $\varepsilon_r=0.95$: hatched = infeasible, white line = $T_4=1000$ K", y=1.01); save(fig, "F01_brayton_surfaces.png")

# F02 P vs A for engineering cases
labs = ["Brayton / B_T", "Brayton / B_C", "", "", "Heat exchanger / Q_max", "Solar / tau-alpha", "Solar / F_R tau-alpha"]
fig, ax = plt.subplots(1, 5, figsize=(14, 3))
for a, c in zip(ax, [1, 2, 5, 6, 7]):
    d = rd(f"st02_sample_{c}.csv"); r = E.iloc[c - 1]
    a.scatter(d["P"], d["A"], s=3, alpha=.4); a.set_xlabel("P"); a.set_ylabel("A = P/B")
    cl = ["R0", "R1", "D2", "I3"][int(r["class_d1"])]
    a.set_title(f"{labs[c-1]}\n{cl}: Rrev(δ)={r['R_rev_delta']:.3f}, Γ={r['Gamma']:.3f}", fontsize=8)
save(fig, "F02_P_vs_A.png")

# F03 rank-reversal map vs x*_P (Carnot benchmark)
M = rd("st09_rrmap.csv"); rr = M["reversal_vs_xP"].values.reshape((nr, nt), order="F")
fig, a = plt.subplots(figsize=(5, 4))
a.contourf(rp, T3, rr, levels=[-.5, .5, 1.5], colors=["#eeeeee", "#d95f02"])
a.contourf(rp, T3, (~feas).astype(float), levels=[.5, 1.5], colors="none", hatches=["///"])
a.plot(*xP, "k*", ms=12); a.set_xlabel("$r_p$"); a.set_ylabel("$T_3$ (K)")
a.set_title(r"Orange: $P<P(x^*_P)$ but $A_T>A_T(x^*_P)$ (rank reversal)" "\n" r"$\varepsilon_r=0.95$", fontsize=8); save(fig, "F03_rank_reversal_map.png")

# F04 Model-2 fronts max[P,A]
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for a, c, t in zip(ax, [1, 2, 5], ["Brayton, Carnot at $T_3$", "Brayton, $B_C$", "Heat exchanger"]):
    d = rd(f"st02_front_{c}.csv"); s = rd(f"st02_sample_{c}.csv"); r = E.iloc[c - 1]
    a.scatter(s["P"], s["A"], s=2, c="0.8", label="feasible sample"); a.scatter(d["P"], d["A"], s=10, c="C3", label="Pareto set, max[P,A]")
    a.plot(r["P_at_xP"], r["A_at_xP"], "k*", ms=11, label="Model 1 optimum"); a.plot(r["P_at_xA"], r["A_at_xA"], "bs", ms=6, label="max A")
    a.set_title(t); a.set_xlabel("P"); a.set_ylabel("A")
    if c != 5:
        a.set_xlim(r["P_at_xA"] - 0.03, r["P_at_xP"] + 0.005); a.set_ylim(r["A_at_xP"] - 0.03, r["A_at_xA"] + 0.005)
ax[0].legend(fontsize=7); save(fig, "F04_model2_fronts.png")

# F05 Model-3 conventional vs benchmark-aware Pareto sets (projected on P-A)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for a, (f, pn, an, t) in zip(ax, [("st03_fronts_brayton_AT.csv", "P", "AT", "Brayton + $A_T$"), ("st03_fronts_hx_A.csv", "P", "A", "Heat exchanger + A"), ("st03_fronts_solar_A2.csv", "P", "A2", "Solar + $A_2$")]):
    d = rd(f); b = d[d.inB == 1]; cc = d[d.inC == 1]
    a.scatter(b[pn], b[an], s=4, c="C1", label="benchmark-aware only"); a.scatter(cc[pn], cc[an], s=4, c="C0", label="conventional (also BA)")
    a.set_title(t); a.set_xlabel("P"); a.set_ylabel(an)
ax[0].legend(fontsize=7); save(fig, "F05_model3_pareto_comparison.png")

# F06 HV convergence
if os.path.exists(os.path.join(R, "st05_hvtrace_M2.csv")):
    T = rd("st05_hvtrace_M2.csv"); fig, a = plt.subplots(figsize=(5, 3.4))
    for c in T.columns: a.plot(np.arange(1, 101) * 100 + 100, T[c], label=c)
    a.set_xlabel("function evaluations"); a.set_ylabel("normalised HV (mean of 30 runs)"); a.set_title("Brayton max[P, $A_C$]"); a.legend(); save(fig, "F06_hv_convergence.png")

# F07 surrogate parity and residuals
Pq = rd("st09_parity.csv"); ms = ["RSM", "GPR", "ANN", "RF", "GBM", "PCGPR"]
fig, ax = plt.subplots(2, 6, figsize=(15, 5))
for j, m in enumerate(ms):
    a = ax[0, j]; a.scatter(Pq["AA"], Pq[f"A_{m}"], s=2); lo, hi = Pq["AA"].min(), Pq["AA"].max(); a.plot([lo, hi], [lo, hi], "k--", lw=.8)
    a.set_title(m.replace("PCGPR", "PC-GPR")); a.set_xlabel("$A_C$ reference model"); a.set_ylabel("$A_C$ surrogate")
    ax[1, j].hist(Pq[f"A_{m}"] - Pq["AA"], bins=40); ax[1, j].set_xlabel("residual in $A_C$")
save(fig, "F07_surrogate_parity_residuals.png")

# F08 Sobol total-order
sb = rd("st06_sobol_brayton.csv"); names = ["P", "$B_T$", "$B_C$", "$A_T$", "$A_C$", "$G_C$", "U", "L"]
fig, ax = plt.subplots(1, 2, figsize=(12, 3.6)); w = .26; x = np.arange(8)
for k, (v, lab) in enumerate(zip(["rp", "T3", "er"], ["$r_p$", "$T_3$", r"$\varepsilon_r$"])):
    y = sb[f"ST_{v}"]; e = [y - sb[f"STlo_{v}"], sb[f"SThi_{v}"] - y]
    ax[0].bar(x + (k - 1) * w, y, w, yerr=e, capsize=2, label=lab)
ax[0].set_xticks(x); ax[0].set_xticklabels(names); ax[0].set_ylabel("total-order Sobol $S_T$"); ax[0].set_title("Brayton (full design box, 95% bootstrap CI)"); ax[0].legend()
sh = rd("st06_sobol_hx.csv"); x = np.arange(5)
for k, (v, lab) in enumerate([("mc", "$\\dot m_c$"), ("A", "$A_{hx}$")]):
    y = sh[f"ST_{v}"]; ax[1].bar(x + (k - .5) * .35, y, .35, yerr=[y - sh[f"STlo_{v}"], sh[f"SThi_{v}"] - y], capsize=2, label=lab)
ax[1].set_xticks(x); ax[1].set_xticklabels(["P=Q", "B=$Q_{max}$", r"A=$\varepsilon$", "U", "L"]); ax[1].set_title("Heat exchanger"); ax[1].legend()
save(fig, "F08_sobol.png")

# F09 gap decomposition
gp = rd("st09_gap_brayton.csv"); parts = list(gp.columns[4:])
fig, a = plt.subplots(figsize=(7, 3.8)); bottom = gp["P"].values.copy()
a.bar(range(3), gp["P"], color="0.3", label=r"$\eta$ achieved")
for k, pcol in enumerate(parts):
    a.bar(range(3), gp[pcol], bottom=bottom, label=pcol); bottom += gp[pcol].values
a.plot(range(3), gp["BT"], "k_", ms=40, mew=2, label="$B_T$ (Carnot at $T_3$)")
a.set_xticks(range(3)); a.set_xticklabels(["max P", "max $A_T$", "max $A_C$"]); a.set_ylabel("efficiency points (fraction of $q_{in}$)")
a.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left"); a.set_title("Exact benchmark-gap decomposition (Gouy-Stodola)"); save(fig, "F09_gap_decomposition.png")

# F10 robustness
sm = rd("st07_samples.csv"); fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
ax[0].hist(sm["P_maxP"], 40, alpha=.6, label="max P design"); ax[0].hist(sm["P_maxAT"], 40, alpha=.6, label="max $A_T$ design"); ax[0].set_xlabel("P under uncertainty"); ax[0].legend(fontsize=7)
ax[1].hist(sm["AT_maxP"], 40, alpha=.6, label="max P design"); ax[1].hist(sm["AT_maxAT"], 40, alpha=.6, label="max $A_T$ design"); ax[1].set_xlabel("$A_T$ under uncertainty"); ax[1].legend(fontsize=7)
d = sm["AT_maxAT"] - sm["AT_maxP"]; ax[2].hist(d, 40); ax[2].axvline(0, c="k"); ax[2].set_xlabel(r"paired $\Delta A_T$ (max-$A_T$ minus max-P)")
ax[2].set_title(f"P(ΔA_T>0) = {np.mean(d > 0):.3f}", fontsize=8); save(fig, "F10_robustness.png")

# F11 optimum displacement in decision space
ro = rd("st07_robust_optima.csv"); fig, a = plt.subplots(figsize=(5, 4))
a.contour(rp, T3, g("T4"), levels=[1000], colors="k", linewidths=1)
for (x, y, m, l) in [(*xP, "r*", "max P"), (*xAT, "bs", "max $A_T$"), (*xAA, "c^", "max $A_C$"),
                     (ro.rp[0], ro.T3[0], "ro", "robust max E[P]"), (ro.rp[1], ro.T3[1], "bo", "robust max E[$A_T$]")]:
    a.plot(x, y, m, ms=9, mfc="none" if "o" in m else None, label=l)
a.annotate("", xy=xAT, xytext=xP, arrowprops=dict(arrowstyle="->")); a.set_xlabel("$r_p$"); a.set_ylabel("$T_3$ (K)")
a.set_title(r"Optimum displacement ($\varepsilon_r=0.95$ at all optima)", fontsize=9); a.legend(fontsize=7); save(fig, "F11_optimum_displacement.png")

# F12 cross-model comparison
sh3 = rd("st03_shift.csv"); cases = [0, 1, 2, 3, 4, 5, 6]; lab = ["Bray/$B_T$", "Bray/$B_C$", "Bray fixed/$B_T$", "Bray fixed/$B_C$", "HX", "Solar/τα", "Solar/$F_R$τα"]
fig, ax = plt.subplots(1, 3, figsize=(13, 3.2))
ax[0].bar(lab, E["R_rev_delta"]); ax[0].set_title(r"$R_{rev}(\delta=1\%)$")
ax[1].bar(lab, E["Gamma"]); ax[1].axhline(.05, c="r", ls="--"); ax[1].set_title(r"$\Gamma$ (dashed = $\gamma_0$)")
ax[2].bar(lab, E["OSI"]); ax[2].set_title("OSI")
for a in ax: a.tick_params(axis="x", rotation=60)
save(fig, "F12_cross_model.png")

# F13 synthetic suite
sy = rd("st01_synthetic.csv"); fig, ax = plt.subplots(1, 3, figsize=(12, 3))
L = ["S1", "S2", "S3", "S4", "S5", "S6"]
ax[0].bar(L, sy["R_rev_delta"]); ax[0].set_title(r"$R_{rev}(\delta)$"); ax[1].bar(L, sy["Gamma"]); ax[1].axhline(.05, c="r", ls="--"); ax[1].set_title(r"$\Gamma$")
ax[2].bar(L, sy["dA_star"], label=r"$\Delta A^*$"); ax[2].bar(L, sy["dP_star"], alpha=.6, label=r"$\Delta P^*$"); ax[2].legend(); ax[2].set_title("optimum regrets")
save(fig, "F13_synthetic_suite.png")
print(sorted(os.listdir(F)))
