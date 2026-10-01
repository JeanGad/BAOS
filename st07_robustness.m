function Rb = st07_robustness(cfg, E)
% ST07_ROBUSTNESS  Monte Carlo uncertainty propagation (m13) at the conventional
% and benchmark optima of the Brayton case; delta-method check for A = P/B;
% chance-constraint violation; robust NSGA-II with common random numbers.
p = m01_define_problem('brayton');
td = struct('field',{'etac','etat','T0','dph'},'type',{'n','n','n','u'},'a',{0.85,0.88,298.15,0.02},'b',{0.01,0.01,5,0.04});
D = [E.xP{1}; E.xA{1}; E.xA{2}];                      % max P, max A_T, max A_A
R = m13_robustness_analysis(p, D, td, cfg.nMC, cfg.seed);
rows = [];
for d = 1:3
  Y = squeeze(R.Y(:,:,d)); P = Y(:,1); BT = Y(:,2); BA = Y(:,3); AT = Y(:,4); AA = Y(:,5);
  dm = @(P,B) sqrt(var(P)/mean(B)^2 + mean(P)^2*var(B)/mean(B)^4 - 2*mean(P)*cov1(P,B)/mean(B)^3);
  % constraint violation probability (T4 > T4max or w_net < w_min)
  viol = 0; for s = 1:cfg.nMC
    t = p.th; for j = 1:numel(td), t.(td(j).field) = R.samples(s,j); end
    o = model_brayton(D(d,:), t); viol = viol + (o.cv > 0);
  end
  rows = [rows; d mean(P) std(P) mean(AT) std(AT) dm(P,BT) mean(AA) std(AA) dm(P,BA) corr2(P,BA) viol/cfg.nMC];
end
write_csv(fullfile(cfg.out, 'st07_mc.csv'), {'design(1=maxP,2=maxAT,3=maxAA)','E_P','sd_P','E_AT','sd_AT','sd_AT_delta','E_AA','sd_AA','sd_AA_delta','corr_P_BA','P_violation'}, rows);
Rb.mc = rows;
% robust formulation with CRN: min[-E(P), -E(A_T), sd(P)] s.t. Pr(violation) <= 5%
rng(cfg.seed + 7); K = 48; Ts = zeros(K, numel(td));
for j = 1:numel(td)
  if td(j).type == 'n', Ts(:,j) = td(j).a + td(j).b*randn(K,1); else, Ts(:,j) = td(j).a + (td(j).b - td(j).a)*rand(K,1); end
end
fun = @(X) robobj(X, p, td, Ts);
res = nsga2(fun, p.lb, p.ub, struct('pop', 80, 'gen', 80, 'seed', cfg.seed));
[~, iP] = min(res.F(:,1)); [~, iA] = min(res.F(:,2));
Rb.robust = [res.X(iP,:) -res.F(iP,1:2) res.F(iP,3); res.X(iA,:) -res.F(iA,1:2) res.F(iA,3)];
write_csv(fullfile(cfg.out, 'st07_robust_front.csv'), {'rp','T3','er','E_P','E_AT','sd_P'}, [res.X -res.F(:,1:2) res.F(:,3)]);
write_csv(fullfile(cfg.out, 'st07_robust_optima.csv'), {'rp','T3','er','E_P','E_AT','sd_P'}, Rb.robust);
% hist data for figure
write_csv(fullfile(cfg.out, 'st07_samples.csv'), {'P_maxP','AT_maxP','P_maxAT','AT_maxAT','P_maxAA','AA_maxAA'}, ...
   [R.Y(:,1,1) R.Y(:,4,1) R.Y(:,1,2) R.Y(:,4,2) R.Y(:,1,3) R.Y(:,5,3)]);
end
function c = cov1(a, b)
c = mean((a - mean(a)).*(b - mean(b)));
end
function r = corr2(a, b)
r = cov1(a,b)/(std(a,1)*std(b,1));
end
function [F, cv] = robobj(X, p, td, Ts)
n = size(X,1); K = size(Ts,1); P = zeros(n,K); AT = P; V = P;
for k = 1:K
  t = p.th; for j = 1:numel(td), t.(td(j).field) = Ts(k,j); end
  o = model_brayton(X, t); P(:,k) = o.P; AT(:,k) = o.AT; V(:,k) = o.cv > 0;
end
F = [-mean(P,2) -mean(AT,2) std(P,0,2)]; cv = max(mean(V,2) - 0.05, 0);
end
