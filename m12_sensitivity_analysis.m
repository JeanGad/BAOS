function S = m12_sensitivity_analysis(fun, lb, ub, N, seed, nboot)
% M12_SENSITIVITY_ANALYSIS  Global Sobol indices (Saltelli 2010 first-order,
% Jansen total-order) with bootstrap 95% intervals. fun: X -> Y (n x q).
d = numel(lb); rng(seed); A = lb + rand(N,d).*(ub - lb); B = lb + rand(N,d).*(ub - lb);
YA = fun(A); YB = fun(B); q = size(YA,2); YAB = zeros(N,q,d);
for i = 1:d, ABi = A; ABi(:,i) = B(:,i); YAB(:,:,i) = fun(ABi); end
[S.S1, S.ST] = est(YA, YB, YAB, 1:N);
S1b = zeros(nboot,q,d); STb = S1b;
for b = 1:nboot, k = randi(N, N, 1); [S1b(b,:,:), STb(b,:,:)] = est(YA, YB, YAB, k); end
S.S1lo = squeeze(prctile1(S1b, 2.5)); S.S1hi = squeeze(prctile1(S1b, 97.5));
S.STlo = squeeze(prctile1(STb, 2.5)); S.SThi = squeeze(prctile1(STb, 97.5));
S.S1lo = reshape(S.S1lo, q, d); S.S1hi = reshape(S.S1hi, q, d); S.STlo = reshape(S.STlo, q, d); S.SThi = reshape(S.SThi, q, d);
end
function [S1, ST] = est(YA, YB, YAB, k)
q = size(YA,2); d = size(YAB,3); S1 = zeros(q,d); ST = zeros(q,d);
V = var([YA(k,:); YB(k,:)], 0, 1);
for i = 1:d
  S1(:,i) = (mean(YB(k,:).*(YAB(k,:,i) - YA(k,:)), 1)./V)';
  ST(:,i) = (0.5*mean((YA(k,:) - YAB(k,:,i)).^2, 1)./V)';
end
end
function p = prctile1(X, q)
X = sort(X, 1); n = size(X,1); i = max(1, min(n, round(q/100*n))); p = X(i,:,:);
end
