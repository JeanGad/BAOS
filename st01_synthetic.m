function S = st01_synthetic(cfg)
% ST01_SYNTHETIC  Level-1/2 validation: BOIT and Model-1 vs Model-2 optimum
% comparison on six problems with analytically known answers.
rng(cfg.seed); n = cfg.nBOIT; X = rand(n, 2);
g = linspace(0,1,201); [a, b] = ndgrid(g, g); Xg = [a(:) b(:)];
S.rows = []; S.names = {};
for c = 1:6
  [P, B, info] = synthetic_suite(X, c);
  r = m04_benchmark_independence(P, B, struct('delta', cfg.delta, 'gamma0', cfg.gamma0), X);
  [Pg, Bg] = synthetic_suite(Xg, c); Ag = Pg./Bg;
  mk = nd_filter([-Pg -Ag]);                          % exact discrete Pareto set of max[P,A]
  [~, iP] = max(Pg + 1e-12*Ag); [~, iA] = max(Ag + 1e-12*Pg);
  OSI = norm(Xg(iA,:) - Xg(iP,:))/sqrt(2);
  % class stability over threshold choices
  cl = {};
  for dl = [0.005 0.01 0.02], for g0 = [0.02 0.05 0.1]
    rr = m04_benchmark_independence(P, B, struct('delta', dl, 'gamma0', g0)); cl{end+1} = rr.class; end, end
  stab = mean(strcmp(cl, r.class));
  % Theorem 2 identity check on 1e5 random pairs
  i1 = randi(n, 1e5, 1); i2 = randi(n, 1e5, 1); A = P./B;
  lhs = P(i1) > P(i2) & A(i1) < A(i2); rhs = P(i1) > P(i2) & B(i1)./B(i2) > P(i1)./P(i2);
  S.names{end+1} = info.name; S.expect{c} = info.expect; S.got{c} = r.class; S.r{c} = r;
  S.rows(end+1,:) = [c, r.CV_B, r.BRI, r.Rrev, r.Rrev_delta, r.tau, r.rhoS, r.Gamma, r.NMI, ...
      r.dA_star, r.dP_star, r.shift, info.shift, nnz(mk), OSI, stab, mean(lhs == rhs), strcmp(r.class, info.expect)];
end
S.header = {'case','CV_B','BRI','R_rev','R_rev_delta','tau','rho_S','Gamma','NMI','dA_star','dP_star', ...
            'shift_detected','shift_expected','n_Pareto_PA','OSI','class_stability','thm2_agreement','class_correct'};
write_csv(fullfile(cfg.out, 'st01_synthetic.csv'), S.header, S.rows);
% store surfaces for figures (S4 and S1)
for c = [1 4]
  [Pg, Bg] = synthetic_suite(Xg, c); write_csv(fullfile(cfg.out, sprintf('st01_surface_S%d.csv', c)), {'x1','x2','P','B','A'}, [Xg Pg Bg Pg./Bg]);
end
end
