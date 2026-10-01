function out = m11_pareto_analysis(F, Fref, ideal, nadir)
% M11_PARETO_ANALYSIS  HV (exact for M<=3, Monte Carlo for M>3), IGD, GD and
% Schott spacing, on objectives normalised by the reference ideal/nadir (ref pt 1.1).
if nargin < 3, ideal = min(Fref,[],1); nadir = max(Fref,[],1); end
sc = max(nadir - ideal, 1e-12); Fn = (F - ideal)./sc; Rn = (Fref - ideal)./sc;
out.HV = hv(Fn, 1.1*ones(1,size(F,2)));
out.IGD = mean(mind(Rn, Fn)); out.GD = mean(mind(Fn, Rn));
if size(Fn,1) > 1
  dd = zeros(size(Fn,1),1);
  for i = 1:size(Fn,1), t = sum(abs(Fn - Fn(i,:)),2); t(i) = inf; dd(i) = min(t); end
  out.SP = std(dd);
else, out.SP = 0; end
end
function d = mind(A, B)
d = zeros(size(A,1),1);
for i = 1:size(A,1), d(i) = sqrt(min(sum((B - A(i,:)).^2, 2))); end
end
