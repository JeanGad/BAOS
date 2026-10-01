function R = st04_surrogates(cfg)
% ST04_SURROGATES  Surrogate comparison for the Brayton case: accuracy (test and
% 5-fold Q2), extrapolation, physical consistency, near-Pareto accuracy, cost;
% then surrogate-assisted MOO with high-fidelity re-evaluation (m15).
p = m01_define_problem('brayton'); lb = p.lb; ub = p.ub; resp = {'P','BA','U','L'};
X = m07_doe(cfg.nTrain, lb, ub, 'lhs', cfg.seed); o = p.eval(X);
rng(cfg.seed + 1); Xt = lb + rand(2000,3).*(ub - lb); ot = p.eval(Xt);
inner = lb + 0.1*(ub - lb); innerU = ub - 0.1*(ub - lb);
Xi = m07_doe(cfg.nTrain, inner, innerU, 'lhs', cfg.seed); oi = p.eval(Xi);
Xo = lb + rand(20000,3).*(ub - lb); Xo = Xo(any(Xo < inner | Xo > innerU, 2),:); Xo = Xo(1:2000,:); oo = p.eval(Xo);
Xc = lb + rand(20000,3).*(ub - lb); oc = p.eval(Xc);
F2 = csvread(fullfile(cfg.out, 'st02_front_2.csv'), 1, 0); Xp = F2(:,1:3); op = p.eval(Xp);    % HF Pareto set of max[P,A_A]
types = {'RSM','GPR','ANN','RF','GBM'}; R.rows = []; R.lab = {};
k = 5; fold = mod(randperm(cfg.nTrain), k) + 1;
for t = 1:numel(types)
  ty = types{t}; printf('  st04: %s\n', ty); fflush(stdout);
  M = struct(); tf = 0; Yt = zeros(2000,4); Yo = Yt; Yc = zeros(size(Xc,1),4); Yp = zeros(size(Xp,1),4); Ycv = zeros(cfg.nTrain,4);
  for r = 1:4
    tic; M.(resp{r}) = m08_surrogate_models('fit', ty, X, o.(resp{r}), lb, ub, cfg.seed); tf = tf + toc;
    tic; Yt(:,r) = m08_surrogate_models('predict', M.(resp{r}), Xt); tp = toc;
    Yc(:,r) = m08_surrogate_models('predict', M.(resp{r}), Xc); Yp(:,r) = m08_surrogate_models('predict', M.(resp{r}), Xp);
    mi = m08_surrogate_models('fit', ty, Xi, oi.(resp{r}), inner, innerU, cfg.seed); Yo(:,r) = m08_surrogate_models('predict', mi, Xo);
    for f = 1:k
      mf = m08_surrogate_models('fit', ty, X(fold ~= f,:), o.(resp{r})(fold ~= f), lb, ub, cfg.seed);
      Ycv(fold == f, r) = m08_surrogate_models('predict', mf, X(fold == f,:));
    end
  end
  R.rows = [R.rows; rowsfor(t, Yt, Ycv, Yo, Yc, Yp, ot, o, oo, oc, op, tf, tp)]; R.lab{t} = ty; R.models{t} = M;
end
% physics-constrained GPR
printf('  st04: PC-GPR\n'); fflush(stdout);
tic; s = pc_gpr('fit', X, o.P, o.BA, o.L, lb, ub, cfg.seed); tf = toc;
[Pp, Bp, Lp] = pc_gpr('predict', s, Xt); tp = 0; Yt = [Pp Bp zeros(2000,1) Lp];
[a1, a2, a3] = pc_gpr('predict', s, Xc); Yc = [a1 a2 zeros(size(a1)) a3]; [a1, a2, a3] = pc_gpr('predict', s, Xp); Yp = [a1 a2 zeros(size(a1)) a3];
si = pc_gpr('fit', Xi, oi.P, oi.BA, oi.L, inner, innerU, cfg.seed); [a1, a2, a3] = pc_gpr('predict', si, Xo); Yo = [a1 a2 zeros(size(a1)) a3];
Ycv = zeros(cfg.nTrain,4);
for f = 1:k
  sf = pc_gpr('fit', X(fold ~= f,:), o.P(fold ~= f), o.BA(fold ~= f), o.L(fold ~= f), lb, ub, cfg.seed);
  [a1, a2, a3] = pc_gpr('predict', sf, X(fold == f,:)); Ycv(fold == f,:) = [a1 a2 zeros(size(a1)) a3];
end
Yt(:,3) = m08_surrogate_models('predict', R.models{2}.U, Xt);   % U from plain GPR
R.rows = [R.rows; rowsfor(6, Yt, Ycv, Yo, Yc, Yp, ot, o, oo, oc, op, tf, tp)]; R.lab{6} = 'PC-GPR'; R.pc = s;
R.head = {'model','resp','R2_test','RMSE_test','MAE_test','MAPE_test','Q2_cv','R2_extrap','RMSE_extrap','MAPE_nearPareto','viol_A_gt_1','viol_L_neg','viol_BA_gt_BT','fit_s','pred_s'};
write_csv(fullfile(cfg.out, 'st04_surrogates.csv'), R.head, R.rows);
save('-v7', fullfile(cfg.out, 'st04_models.mat'), 'R', 'X', 'o');
end
function rows = rowsfor(t, Yt, Ycv, Yo, Yc, Yp, ot, o, oo, oc, op, tf, tp)
names = {'P','BA','U','L'}; rows = [];
for r = [1 2 3 4 5]
  if r <= 4
    y = ot.(names{r}); yh = Yt(:,r); yc = o.(names{r}); ycv = Ycv(:,r); yo = oo.(names{r}); yho = Yo(:,r); ypf = op.(names{r}); yph = Yp(:,r);
  else   % derived attainment A = P_hat / B_hat
    y = ot.AA; yh = Yt(:,1)./Yt(:,2); yc = o.AA; ycv = Ycv(:,1)./Ycv(:,2); yo = oo.AA; yho = Yo(:,1)./Yo(:,2); ypf = op.AA; yph = Yp(:,1)./Yp(:,2);
  end
  if t == 6 && r == 3, rows = [rows; t r nan(1,13)]; continue; end
  mt = m09_surrogate_validation(y, yh); mc = m09_surrogate_validation(yc, ycv); me = m09_surrogate_validation(yo, yho); mp = m09_surrogate_validation(ypf, yph);
  vA = mean(Yc(:,1) > Yc(:,2)); vL = mean(Yc(:,4) < 0); vB = mean(Yc(:,2) > oc.BT);
  rows = [rows; t r mt.R2 mt.RMSE mt.MAE mt.MAPE mc.R2 me.R2 me.RMSE mp.MAPE vA vL vB tf tp];
end
end
