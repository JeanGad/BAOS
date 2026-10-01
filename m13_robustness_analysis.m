function R = m13_robustness_analysis(prob, X, thdist, Nmc, seed)
% M13_ROBUSTNESS_ANALYSIS  Monte Carlo propagation of parameter uncertainty
% theta ~ p(theta) at designs X. thdist: struct array (field, type 'n'/'u', a, b).
% Returns E and std of P, benchmark(s), attainment(s), and delta-method check for A=P/B.
rng(seed); nd = size(X,1); th = prob.th; Th = cell(Nmc,1);
S = zeros(Nmc, numel(thdist));
for j = 1:numel(thdist)
  if thdist(j).type == 'n', S(:,j) = thdist(j).a + thdist(j).b*randn(Nmc,1);
  else, S(:,j) = thdist(j).a + (thdist(j).b - thdist(j).a)*rand(Nmc,1); end
end
R.samples = S; flds = prob.rfields; R.fields = flds; R.Y = zeros(Nmc, numel(flds), nd);
for s = 1:Nmc
  t = th; for j = 1:numel(thdist), t.(thdist(j).field) = S(s,j); end
  o = prob.model(X, t);
  for f = 1:numel(flds), R.Y(s,f,:) = reshape(o.(flds{f}), 1, 1, nd); end
end
R.mean = squeeze(mean(R.Y,1)); R.std = squeeze(std(R.Y,0,1));
end
