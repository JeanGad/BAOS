import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9.5, "savefig.dpi": 300})
# Fig 2: Theorem 1 (gradient alignment) illustration with labelled arrows
x1, x2 = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 1, 401))
P = 1 - (x1 - .7)**2 - (x2 - .7)**2; B = 1.5 + x1; A = P/B
Pf, Af, X1, X2 = P.ravel(), A.ravel(), x1.ravel(), x2.ravel()
o = np.argsort(-Pf); best = -np.inf; m = np.zeros(Pf.size, bool)
for i in o:
    if Af[i] > best + 1e-12: m[i] = True; best = Af[i]
fig, ax = plt.subplots(figsize=(5, 4.4))
ax.contour(x1, x2, P, 10, colors="#0072B2", linewidths=.7); ax.contour(x1, x2, A, 10, colors="#D55E00", linestyles="--", linewidths=.7)
ax.plot(X1[m], X2[m], "k-", lw=3, label="Pareto set of max[P, A]")
def arr(x, y, dx, dy, c, t, off):
    ax.annotate("", xy=(x+dx, y+dy), xytext=(x, y), arrowprops=dict(arrowstyle="->", color=c, lw=1.5))
    ax.text(x+dx+off[0], y+dy+off[1], t, color=c, fontsize=8.5)
for xa, lab, pos in [(0.55, "0 ≤ κ ≤ A:\ncritical", (0.50, 0.58)), (0.25, "κ = 0.9 > A = 0.46:\nnot critical", (0.03, 0.58))]:
    gP = -2*(xa-.7)
    arr(xa, .7, 0.25*gP, 0, "#0072B2", "∇P", (0.01, -0.035)); arr(xa, .74, 0.12, 0.0, "#D55E00", "∇B", (0.01, 0.01))
    ax.text(pos[0], pos[1], lab, fontsize=7.5)
xo, yo = 0.45, 0.3; gP = np.array([-2*(xo-.7), -2*(yo-.7)])
arr(xo, yo, 0.25*gP[0], 0.25*gP[1], "#0072B2", "∇P", (0.01, 0.0)); arr(xo, yo, 0.12, 0, "#D55E00", "∇B", (0.01, -0.03))
ax.text(xo-0.05, yo-0.08, "not parallel: not critical", fontsize=7.5)
ax.plot([], [], color="#0072B2", label="contours of P"); ax.plot([], [], color="#D55E00", ls="--", label="contours of A")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$"); ax.legend(fontsize=7.5, loc="lower left")
fig.tight_layout(); fig.savefig("fig02_theorem_illustration.png", bbox_inches="tight"); plt.close(fig)
# Composite figure: three stacked panels, larger
g = np.loadtxt("hs_grid.csv", delimiter=","); k = 201
F = g[:,0].reshape(k,k, order="F"); T = g[:,1].reshape(k,k, order="F")
fig, ax = plt.subplots(3, 1, figsize=(5.6, 11))
for a, j, t in zip(ax, [2, 4, 5], ["(a) $P=K^*/\\rho$ (GPa cm$^3$ g$^{-1}$)", "(b) $A_{HS}=K^*/K^{U}_{HS}$", "(c) bound-span attainment $(K^*-K^{L}_{HS})/(K^{U}_{HS}-K^{L}_{HS})$"]):
    cs = a.contourf(F, T, g[:,j].reshape(k,k, order="F"), 20, cmap="viridis"); fig.colorbar(cs, ax=a)
    a.plot(0.60, 0.50, "r*", ms=14, mec="w", label="max P", clip_on=False); a.set_xlabel("stiff-phase fraction $f$"); a.set_ylabel(r"connectivity index $\theta$"); a.set_title(t, fontsize=10)
ax[1].plot(0.05, 0.50, "ws", ms=9, mec="k", label="max $A_{HS}$", clip_on=False); ax[2].plot(0.60, 0.485, "c^", ms=9, mec="k", label="max bound-span", clip_on=False)
for a in ax: a.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig("fig_composite_v2.png", bbox_inches="tight"); plt.close(fig)
# Robustness figure with three panels
s1 = pd.read_csv("study1.csv"); s2 = pd.read_csv("study2.csv"); s2b = pd.read_csv("study2b.csv"); gr = pd.read_csv("graded.csv")
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
mc = s1.groupby("n").class_acc.mean(); s6 = s1[s1["case"] == 6].set_index("n").shift_acc
ax[0].plot(mc.index, mc.values, "o-", color="#0072B2", label="class correct (S1–S6)"); ax[0].plot(s6.index, s6.values, "s--", color="#D55E00", label="sample-based shift verdict, S6")
ax[0].set_xscale("log"); ax[0].set_xlabel("sample size n"); ax[0].set_ylabel("fraction correct"); ax[0].set_ylim(-.05, 1.05); ax[0].legend(fontsize=7.5); ax[0].set_title("(a) sample size")
x = [100*v for v in sorted(s2.sigma.unique())]
ax[1].plot(x, s2.groupby("sigma").class_acc_d1.mean().values, "o-", color="#0072B2", label="δ = 1 %, r$_{tol}$ = 0")
ax[1].plot(x, s2.groupby("sigma").class_acc_d3s.mean().values, "s-", color="#E69F00", label="δ = max(1 %, 3σ), r$_{tol}$ = 0")
ax[1].plot(x, s2b.groupby("sigma").acc_rtol005.mean().values, "^-", color="#009E73", label="δ = max(1 %, 3σ), r$_{tol}$ = 0.5 %")
ax[1].set_xlabel("relative noise σ (%)"); ax[1].set_ylabel("class correct (S1–S6)"); ax[1].set_ylim(-.05, 1.05); ax[1].legend(fontsize=7.5); ax[1].set_title("(b) evaluation noise")
for (sg, rh), mk, c, lab in [((0, 0), "o-", "#0072B2", "no noise"), ((0.01, 0), "s--", "#E69F00", "σ = 1 %, independent"), ((0.01, 0.8), "^:", "#009E73", "σ = 1 %, correlated (ρ = 0.8)")]:
    d = gr[(gr.sigma == sg) & (gr.rho == rh)]; d = d[d.u_pop > 0]
    ax[2].plot(d.u_pop, d.frac_I3, mk, color=c, label=lab)
ax[2].axvline(0.05, color="k", lw=.7, ls="--"); ax[2].text(0.055, 0.5, "γ₀", fontsize=9)
ax[2].set_xscale("log"); ax[2].set_xlabel("unexplained variance share of A given P"); ax[2].set_ylabel("fraction classified I3"); ax[2].set_ylim(-.05, 1.05); ax[2].legend(fontsize=7.5); ax[2].set_title("(c) graded independence (unseen family)")
fig.tight_layout(); fig.savefig("fig_robust_v2.png", bbox_inches="tight"); plt.close(fig)
print("ok")
