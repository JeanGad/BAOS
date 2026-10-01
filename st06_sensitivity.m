function S = st06_sensitivity(cfg, E)
% ST06_SENSITIVITY  Global Sobol (first/total order, bootstrap 95% CI) for
% P, B_T, B_A, A_T, A_A, G_A, U, L (Brayton) and P, B, A, U, L (heat exchanger);
% local log-elasticities at x*_P and x*_A with a numerical check of
% Proposition 5: e(A) = e(P) - e(B).
p = m01_define_problem('brayton'); nb = {'P','BT','BA','AT','AA','G','U','L'};
f = @(X) brayY(p.eval(X));
S.bray = m12_sensitivity_analysis(f, p.lb, p.ub, cfg.nSobol, cfg.seed, 200);
rows = []; for i = 1:numel(nb), rows = [rows; i S.bray.S1(i,:) S.bray.S1lo(i,:) S.bray.S1hi(i,:) S.bray.ST(i,:) S.bray.STlo(i,:) S.bray.SThi(i,:)]; end
write_csv(fullfile(cfg.out, 'st06_sobol_brayton.csv'), {'resp','S1_rp','S1_T3','S1_er','S1lo_rp','S1lo_T3','S1lo_er','S1hi_rp','S1hi_T3','S1hi_er', ...
   'ST_rp','ST_T3','ST_er','STlo_rp','STlo_T3','STlo_er','SThi_rp','SThi_T3','SThi_er'}, rows);
q = m01_define_problem('hx'); g = @(X) hxY(q.eval(X));
S.hx = m12_sensitivity_analysis(g, q.lb, q.ub, cfg.nSobol, cfg.seed, 200);
rows = []; for i = 1:5, rows = [rows; i S.hx.S1(i,:) S.hx.ST(i,:) S.hx.STlo(i,:) S.hx.SThi(i,:)]; end
write_csv(fullfile(cfg.out, 'st06_sobol_hx.csv'), {'resp','S1_mc','S1_A','ST_mc','ST_A','STlo_mc','STlo_A','SThi_mc','SThi_A'}, rows);
% local elasticities (Proposition 5)
rows = [];
for c = 1:2
  bn = {'BT','BA'}; an = {'AT','AA'};
  for w = 1:2
    if w == 1, x = E.xP{c}; else, x = E.xA{c}; end
    x(3) = min(x(3), 0.95*(1 - 1e-4));                          % keep central difference inside bounds
    h = @(X) pick3(p.eval(X), bn{c}, an{c}); El = local_elasticity(h, x, 1e-5);
    rows = [rows; c w El(1,:) El(2,:) El(3,:) max(abs(El(3,:) - (El(1,:) - El(2,:))))];
  end
end
write_csv(fullfile(cfg.out, 'st06_elasticity_brayton.csv'), {'bench(1=BT,2=BA)','point(1=xP,2=xA)','eP_rp','eP_T3','eP_er','eB_rp','eB_T3','eB_er','eA_rp','eA_T3','eA_er','max_identity_residual'}, rows);
S.elas = rows;
end
function Y = brayY(o)
Y = [o.P o.BT o.BA o.AT o.AA 1 - o.AA o.U o.L];
end
function Y = hxY(o)
Y = [o.P o.B o.A o.U o.L];
end
function y = pick3(o, bn, an)
y = [o.P o.(bn) o.(an)];
end
