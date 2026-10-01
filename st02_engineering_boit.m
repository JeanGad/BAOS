function E = st02_engineering_boit(cfg)
% ST02_ENGINEERING_BOIT  BOIT, exact Model-1 vs Model-2 optima, Pareto set of
% max[P,A] and numerical check of Theorem 4 for each engineering case.
cases = {'brayton','BT',[],'Brayton / Carnot at T3 (B_T)'; 'brayton','BA',[],'Brayton / heat-exergy limit (B_A)'; ...
         'brayton','BT',1400,'Brayton, T3 fixed / B_T'; 'brayton','BA',1400,'Brayton, T3 fixed / B_A'; ...
         'hx','B',[],'Heat exchanger / Q_max'; 'solar','B1',[],'Solar / tau*alpha'; 'solar','B2',[],'Solar / F_R*tau*alpha'};
E.rows = []; E.labels = cases(:,4);
for c = 1:size(cases,1)
  p = m01_define_problem(cases{c,1}); lb = p.lb; ub = p.ub; bn = cases{c,2};
  if ~isempty(cases{c,3}), lb(2) = cases{c,3}; ub(2) = cases{c,3}; end
  ev = @(X) addA(p.eval(X), bn);
  % feasible space-filling sample for BOIT
  rng(cfg.seed + c); X = lb + rand(4*cfg.nBOIT, numel(lb)).*(ub - lb); o = ev(X); X = X(o.cv <= 0,:); X = X(1:cfg.nBOIT,:); o = ev(X);
  r = m04_benchmark_independence(o.P, o.B, struct('delta', cfg.delta, 'gamma0', cfg.gamma0), X);
  r3 = m04_benchmark_independence(o.P, o.B, struct('delta', cfg.delta/10, 'gamma0', cfg.gamma0), X);
  % dense grid
  d = numel(lb); k = ternary(d - nnz(ub == lb) >= 3, 41, 201); Xg = m07_doe(k^d, lb, ub, 'grid', 1);
  if nnz(ub == lb) > 0, Xg = unique(Xg, 'rows'); end
  % exact (lexicographic, tie tolerance 1e-9) and delta-resolved optima
  [xP, oP] = lex_opt(ev, 'P', 'A', lb, ub, Xg, 1e-9); [xA, oA] = lex_opt(ev, 'A', 'P', lb, ub, Xg, 1e-9);
  [~, oPd] = lex_opt(ev, 'P', 'A', lb, ub, Xg, cfg.delta);
  dA = (oA.A - oP.A)/oA.A; dP = (oP.P - oA.P)/oP.P; dAd = (oA.A - oPd.A)/oA.A;
  sc = ub - lb; sc(sc == 0) = 1; fr = ub > lb; OSI = norm((xA(fr) - xP(fr))./sc(fr))/sqrt(nnz(fr));
  % Pareto set of max[P,A] by NSGA-II and Theorem 4 check on interior, inactive-constraint points
  fun = @(Z) wrapPA(ev, Z);
  if nnz(fr) >= 1
    res = nsga2(fun, lb, ub, struct('pop', 100, 'gen', 150, 'seed', cfg.seed));
    [t4n, t4ok, t4cos] = thm4check(p, bn, res.X, lb, ub);
  else, res.X = xP; res.F = [-oP.P -oP.A]; t4n = 0; t4ok = NaN; t4cos = NaN; end
  write_csv(fullfile(cfg.out, sprintf('st02_front_%d.csv', c)), [p.vars(:)' {'P','A'}], [res.X -res.F]);
  og = ev(Xg); fe = og.cv <= 0; mk = false(size(fe)); mk(fe) = nd_filter([-og.P(fe) -og.A(fe)]);
  write_csv(fullfile(cfg.out, sprintf('st02_sample_%d.csv', c)), [p.vars(:)' {'P','B','A'}], [X o.P o.B o.A]);
  E.r{c} = r; E.xP{c} = xP; E.xA{c} = xA; E.oP{c} = oP; E.oA{c} = oA;
  E.rows(end+1,:) = [c, r.validity, r.CV_B, r.BRI, r.Rrev, r.Rrev_delta, r.Rrev_top, r.tau, r.rhoS, r.Gamma, r.Gamma_top, r.NMI, ...
      cls(r.class), cls(r3.class), oP.P, oP.A, oA.P, oA.A, dA, dP, dAd, OSI, nnz(mk), size(res.F,1), t4n, t4ok, t4cos, [xP nan(1,3-numel(xP))], [xA nan(1,3-numel(xA))]];
  E.cls{c} = r.class; E.cls3{c} = r3.class;
end
E.header = [{'case','validity','CV_B','BRI','R_rev','R_rev_delta','R_rev_top','tau','rho_S','Gamma','Gamma_top','NMI', ...
     'class_d1','class_d01','P_at_xP','A_at_xP','P_at_xA','A_at_xA','dA_exact','dP_exact','dA_delta','OSI','nPareto_grid','nPareto_nsga', ...
     'thm4_n','thm4_frac_ok','thm4_median_cos'} strcat('xP_', {'1','2','3'}) strcat('xA_', {'1','2','3'})];
E.rows = pad(E.rows, numel(E.header));
write_csv(fullfile(cfg.out, 'st02_engineering_boit.csv'), E.header, E.rows);
end
function o = addA(o, bn)
o.B = o.(bn); o.A = o.P./o.B;
end
function [F, cv] = wrapPA(ev, Z)
o = ev(Z); F = [-o.P -o.A]; cv = o.cv;
end
function k = cls(s)
k = find(strcmp(s, {'R0','R1','D2','I3'})) - 1;
end
function y = ternary(c, a, b)
if c, y = a; else, y = b; end
end
function M = pad(M, n)
if size(M,2) < n, M = [M nan(size(M,1), n - size(M,2))]; end
end
function [n, fok, mcos] = thm4check(p, bn, X, lb, ub)
% Theorem 4 / Corollary 4.1: at Pareto-critical points of max[P,A], the gradients
% of P and B projected onto the tangent space of the active constraints satisfy
% Pi*grad(P) = kappa * Pi*grad(B) with 0 <= kappa <= A. Bounds at their limits are
% handled by removing those coordinates; active nonlinear constraints by projection.
h = 1e-6; ok = []; cs = [];
for i = 1:size(X,1)
  x = X(i,:); o = p.eval(x); if o.cv > 0, continue; end
  fr = find((x - lb) > 1e-4*(ub - lb) & (ub - x) > 1e-4*(ub - lb)); if isempty(fr), continue; end
  act = []; if ~isempty(o.g), act = find(o.g./p.gscale > -1e-3); end
  nf = numel(fr); gP = zeros(nf,1); gB = gP; Gc = zeros(nf, numel(act));
  for j = 1:nf
    e = zeros(size(x)); e(fr(j)) = h*(ub(fr(j)) - lb(fr(j))); op = p.eval(x + e); om = p.eval(x - e);
    gP(j) = (op.P - om.P)/(2*h); gB(j) = (op.(bn) - om.(bn))/(2*h);
    if ~isempty(act), Gc(j,:) = (op.g(act) - om.g(act))./p.gscale(act)/(2*h); end
  end
  if numel(act) >= nf, continue; end                       % vertex: no tangent space left
  if isempty(act), Pi = eye(nf); else, Pi = eye(nf) - Gc*pinv(Gc); end
  a = Pi*gP; b = Pi*gB; A = o.P/o.(bn);
  kap = (a'*b)/(b'*b + eps); cs(end+1) = (a'*b)/(norm(a)*norm(b) + eps);
  resid = norm(a - kap*b)/(norm(a) + norm(b) + eps);
  ok(end+1) = (kap >= -1e-6) && (kap <= A + 1e-6) && resid < 1e-2;
end
n = numel(ok); if n == 0, fok = NaN; mcos = NaN; else, fok = mean(ok); mcos = median(cs); end
end
