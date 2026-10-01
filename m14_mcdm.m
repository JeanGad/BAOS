function out = m14_mcdm(D, benefit, w)
% M14_MCDM  TOPSIS with vector normalisation. D: alternatives x criteria,
% benefit: logical (true = larger is better), w: weights (default: entropy).
% Also returns rank stability under +/-10% and +/-20% one-at-a-time weight
% perturbations (fraction of perturbations that keep the same top choice).
[n, m] = size(D);
if nargin < 3 || isempty(w)
  Pn = (D - min(D) + eps)./sum(D - min(D) + eps, 1); E = -sum(Pn.*log(Pn), 1)/log(n);
  w = (1 - E)/sum(1 - E);
end
out.w = w(:)'; [out.C, out.best] = topsis(D, benefit, out.w);
keep = []; for f = [0.9 1.1 0.8 1.2], for j = 1:m
  w2 = out.w; w2(j) = w2(j)*f; w2 = w2/sum(w2); [~, b] = topsis(D, benefit, w2); keep(end+1) = (b == out.best); end, end
out.stability = mean(keep);
end
function [C, best] = topsis(D, benefit, w)
V = D./sqrt(sum(D.^2,1)).*w; ip = max(V,[],1); an = min(V,[],1);
ip(~benefit) = min(V(:,~benefit),[],1); an(~benefit) = max(V(:,~benefit),[],1);
dp = sqrt(sum((V - ip).^2,2)); dn = sqrt(sum((V - an).^2,2)); C = dn./(dp + dn + eps); [~, best] = max(C);
end
