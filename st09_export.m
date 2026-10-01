function st09_export(cfg, E)
% ST09_EXPORT  Figure data: Brayton surfaces (P, B, A, G) at eps_r = 0.95,
% rank-reversal map, gap decomposition at key designs, surrogate parity data.
p = m01_define_problem('brayton');
[a, b] = ndgrid(linspace(2, 24, 111), linspace(1100, 1700, 121)); X = [a(:) b(:) 0.95*ones(numel(a),1)]; o = p.eval(X);
write_csv(fullfile(cfg.out, 'st09_surface_brayton.csv'), {'rp','T3','P','BT','BA','AT','AA','U','L','T4','cv'}, [X(:,1:2) o.P o.BT o.BA o.AT o.AA o.U o.L o.T4 o.cv]);
% rank-reversal map relative to the conventional optimum x*_P (Carnot benchmark)
x0 = E.xP{1}; o0 = p.eval(x0);
rr = (o.P < o0.P & o.AT > o0.AT); write_csv(fullfile(cfg.out, 'st09_rrmap.csv'), {'rp','T3','reversal_vs_xP','dominated'}, [X(:,1:2) rr (o.P <= o0.P & o.AT <= o0.AT)]);
% gap decomposition at x*_P, x*_A(B_T), x*_A(B_A)
D = [E.xP{1}; E.xA{1}; E.xA{2}]; od = p.eval(D);
write_csv(fullfile(cfg.out, 'st09_gap_brayton.csv'), [{'design','P','BT','BA'} od.partnames], [(1:3)' od.P od.BT od.BA od.parts]);
q = m01_define_problem('hx'); D = [E.xP{5}; E.xA{5}]; oq = q.eval(D);
exIn = oq.U + oq.parts(:,1);
write_csv(fullfile(cfg.out, 'st09_gap_hx.csv'), {'design','Q','Qmax','eps','ex_gained','exD_HT','W_pump'}, [(1:2)' oq.P oq.B oq.A oq.U oq.parts]);
s = m01_define_problem('solar'); [a2, b2] = ndgrid(linspace(s.lb(1), s.ub(1), 5), [298.15 318.15 338.15]); Ds = [a2(:) b2(:)]; os = s.eval(Ds);
write_csv(fullfile(cfg.out, 'st09_gap_solar.csv'), {'m','Ti','eta','ta','heat_removal_gap','thermal_loss_gap','closure'}, ...
   [Ds os.P os.B1 os.B1 - os.B2 os.B2 - os.P os.B1 - os.P - ((os.B1 - os.B2) + (os.B2 - os.P))]);
% surrogate parity (test set) for the five learners + PC-GPR on P and A_A
L = load(fullfile(cfg.out, 'st04_models.mat')); R = L.R; rng(cfg.seed + 1); Xt = p.lb + rand(2000,3).*(p.ub - p.lb); ot = p.eval(Xt);
M = [ot.P ot.AA];
for t = 1:5
  Pp = m08_surrogate_models('predict', R.models{t}.P, Xt); Bp = m08_surrogate_models('predict', R.models{t}.BA, Xt); M = [M Pp Pp./Bp];
end
[Pp, Bp] = pc_gpr('predict', R.pc, Xt); M = [M Pp Pp./Bp];
write_csv(fullfile(cfg.out, 'st09_parity.csv'), {'P','AA','P_RSM','A_RSM','P_GPR','A_GPR','P_ANN','A_ANN','P_RF','A_RF','P_GBM','A_GBM','P_PCGPR','A_PCGPR'}, M);
end
