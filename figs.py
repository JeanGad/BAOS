import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "savefig.dpi": 300})
# ---- Fig A: illustration of Theorem 4 (interior Pareto path) ----
x1, x2 = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 1, 401))
P = 1 - (x1 - .7)**2 - (x2 - .7)**2; B = 1.5 + x1; A = P/B
Pf, Af, X1, X2 = P.ravel(), A.ravel(), x1.ravel(), x2.ravel()
o = np.argsort(-Pf); best = -np.inf; m = np.zeros(Pf.size, bool)
for i in o:
    if Af[i] > best + 1e-12: m[i] = True; best = Af[i]
fig, ax = plt.subplots(figsize=(4.6, 4))
c1 = ax.contour(x1, x2, P, 10, colors="C0", linewidths=.7); c2 = ax.contour(x1, x2, A, 10, colors="C3", linestyles="--", linewidths=.7)
ax.plot(X1[m], X2[m], "k.", ms=3, label="Pareto set of max[P, A]")
for xa in [0.25, 0.45, 0.65]:
    gP = np.array([-2*(xa-.7), 0.0]); gB = np.array([1.0, 0.0])
    ax.annotate("", xy=(xa+0.25*gP[0], .7+0.25*gP[1]), xytext=(xa, .7), arrowprops=dict(arrowstyle="->", color="C0", lw=1.4))
    ax.annotate("", xy=(xa+0.12*gB[0], .7+0.12*gB[1]+0.03), xytext=(xa, .73), arrowprops=dict(arrowstyle="->", color="C3", lw=1.4))
xo, yo = 0.45, 0.3; gP = np.array([-2*(xo-.7), -2*(yo-.7)])
ax.annotate("", xy=(xo+0.25*gP[0], yo+0.25*gP[1]), xytext=(xo, yo), arrowprops=dict(arrowstyle="->", color="C0", lw=1.4))
ax.annotate("", xy=(xo+0.12, yo), xytext=(xo, yo), arrowprops=dict(arrowstyle="->", color="C3", lw=1.4))
ax.text(xo+0.02, yo-0.07, "off the set:\n∇P not parallel to ∇B", fontsize=7.5)
ax.plot([], [], color="C0", label="contours of P; arrows ∇P"); ax.plot([], [], color="C3", ls="--", label="contours of A; arrows ∇B")
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$"); ax.legend(fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig("figR1_theorem4_illustration.png", bbox_inches="tight"); plt.close(fig)
# ---- Fig B: composite case ----
g = np.loadtxt("hs_grid.csv", delimiter=","); k = 201
F = g[:,0].reshape(k,k, order="F"); T = g[:,1].reshape(k,k, order="F")
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for a, j, t in zip(ax, [2, 4, 5], ["P: specific bulk modulus (GPa cm$^3$/g)", "$A_{HS}=K^*/K^{HS+}$", "bound-span attainment $(K^*-K^{HS-})/(K^{HS+}-K^{HS-})$"]):
    cs = a.contourf(F, T, g[:,j].reshape(k,k, order="F"), 20, cmap="viridis"); fig.colorbar(cs, ax=a)
    a.plot(0.60, 0.50, "r*", ms=12, label="max P"); a.set_xlabel("stiff-phase fraction $f$"); a.set_ylabel(r"connectivity index $\theta$"); a.set_title(t, fontsize=8.5)
ax[1].plot(0.05, 0.50, "ws", ms=7, mec="k", label="max $A_{HS}$"); ax[2].plot(0.60, 0.50, "c^", ms=7, mec="k", label="max bound-span")
for a in ax: a.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig("figR2_composite.png", bbox_inches="tight"); plt.close(fig)
# ---- Fig C: robustness ----
s1 = pd.read_csv("study1.csv"); s2 = pd.read_csv("study2.csv"); s2b = pd.read_csv("study2b.csv")
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
mean_cls = s1.groupby("n").class_acc.mean(); s6 = s1[s1["case"] == 6].set_index("n").shift_acc
ax[0].plot(mean_cls.index, mean_cls.values, "o-", label="class correct (mean of S1–S6)")
ax[0].plot(s6.index, s6.values, "s--", label="sample-based shift verdict, S6 (corner optima)")
ax[0].axhline(1, color="k", lw=.5); ax[0].set_xscale("log"); ax[0].set_xlabel("feasible sample size n"); ax[0].set_ylabel("fraction correct (20 replicates)"); ax[0].legend(fontsize=7); ax[0].set_ylim(-.05, 1.05)
sg = sorted(s2.sigma.unique()); x = [100*v for v in sg]
ax[1].plot(x, s2.groupby("sigma").class_acc_d1.mean().values, "o-", label="δ = 1 %, zero tolerance")
ax[1].plot(x, s2.groupby("sigma").class_acc_d3s.mean().values, "s-", label="δ = max(1 %, 3σ), zero tolerance")
ax[1].plot(x, s2b.groupby("sigma").acc_rtol005.mean().values, "^-", label="δ = max(1 %, 3σ), r_tol = 0.5 %")
ax[1].set_xlabel("relative evaluation noise σ (%)"); ax[1].set_ylabel("class correct (mean of S1–S6)"); ax[1].legend(fontsize=7); ax[1].set_ylim(-.05, 1.05)
fig.tight_layout(); fig.savefig("figR3_boit_robustness.png", bbox_inches="tight"); plt.close(fig)
print("ok")
