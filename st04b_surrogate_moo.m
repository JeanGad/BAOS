function V = st04b_surrogate_moo(cfg, E)
% ST04B_SURROGATE_MOO  Surrogate-assisted NSGA-II for max[P, A_A] and max[P, A_T]
% (GPR surrogates for P, B_A, T4, w_net; B_T is analytic), followed by
% high-fidelity re-evaluation of every Pareto candidate (m15) and GPR
% uncertainty propagation into A = P/B (delta method) with coverage check.
p = m01_define_problem('brayton'); lb = p.lb; ub = p.ub;
X = m07_doe(cfg.nTrain, lb, ub, 'lhs', cfg.seed); o = p.eval(X);
S.P = m08_surrogate_models('fit', cfg.surrogate, X, o.P, lb, ub, cfg.seed);
S.BA = m08_surrogate_models('fit', cfg.surrogate, X, o.BA, lb, ub, cfg.seed);
S.T4 = m08_surrogate_models('fit', cfg.surrogate, X, o.T4, lb, ub, cfg.seed);
S.U = m08_surrogate_models('fit', cfg.surrogate, X, o.U, lb, ub, cfg.seed);
rows = []; bn = {'AA','AT'};
for c = 1:2
  fun = @(Z) surfun(Z, S, p, bn{c});
  res = nsga2(fun, lb, ub, struct('pop', 100, 'gen', 150, 'seed', cfg.seed));
  oh = p.eval(res.X); Ps = -res.F(:,1); As = -res.F(:,2);
  if c == 1, Bs = m08_surrogate_models('predict', S.BA, res.X); Bh = oh.BA; Ah = oh.AA; else, Bs = 1 - p.th.T0./res.X(:,2); Bh = oh.BT; Ah = oh.AT; end
  eP = 100*abs(oh.P - Ps)./oh.P; eB = 100*abs(Bh - Bs)./Bh; eA = 100*abs(Ah - As)./Ah;
  % HF front of the same problem (from st02) for IGD in HF objective space
  Fhf = csvread(fullfile(cfg.out, sprintf('st02_front_%d.csv', 3 - c)), 1, 0); Fhf = -Fhf(:,4:5);
  feas = oh.cv <= 0; Fs_hf = [-oh.P(feas) -Ah(feas)];
  q = m11_pareto_analysis(Fs_hf(nd_filter(Fs_hf),:), Fhf);
  rows = [rows; c size(res.X,1) mean(eP) max(eP) mean(eB) max(eB) mean(eA) max(eA) mean(~feas) q.IGD q.GD];
  write_csv(fullfile(cfg.out, sprintf('st04b_surfront_%s.csv', bn{c})), {'rp','T3','er','P_sur','A_sur','P_HF','A_HF','B_sur','B_HF','cv_HF'}, [res.X Ps As oh.P Ah Bs Bh oh.cv]);
end
write_csv(fullfile(cfg.out, 'st04b_hf_validation.csv'), {'case(1=AA,2=AT)','nPareto','eP_mean','eP_max','eB_mean','eB_max','eA_mean','eA_max','frac_HF_infeasible','IGD_vs_HF','GD_vs_HF'}, rows);
V.rows = rows;
% uncertainty propagation for A_A on an independent test set
rng(cfg.seed + 3); Xt = lb + rand(2000,3).*(ub - lb); ot = p.eval(Xt);
[mP, sP] = m08_surrogate_models('predict', S.P, Xt); [mB, sB] = m08_surrogate_models('predict', S.BA, Xt);
mA = mP./mB; sA = sqrt(sP.^2./mB.^2 + mP.^2.*sB.^2./mB.^4);          % independent GPs: Cov(P,B) = 0
z = abs(ot.AA - mA)./max(sA, eps);
V.cov = [mean(z <= 1.96) median(sA) median(abs(ot.AA - mA))];
write_csv(fullfile(cfg.out, 'st04b_uq.csv'), {'coverage95','median_sd_A','median_abs_err_A'}, V.cov);
end
function [F, cv] = surfun(Z, S, p, bn)
P = m08_surrogate_models('predict', S.P, Z);
if strcmp(bn, 'AA'), B = m08_surrogate_models('predict', S.BA, Z); else, B = 1 - p.th.T0./Z(:,2); end
T4 = m08_surrogate_models('predict', S.T4, Z); U = m08_surrogate_models('predict', S.U, Z);
F = [-P -P./B]; cv = max(T4 - p.th.T4max, 0)/p.th.T4max + max(p.th.wmin - U, 0)/p.th.wmin;
end
