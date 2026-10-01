function r = m04_benchmark_independence(P, B, opt, X)
% M04_BENCHMARK_INDEPENDENCE  Benchmark Objective Independence Test (BOIT).
% P, B : performance and benchmark over a space-filling sample of the FEASIBLE set
% opt  : delta (relative resolution), gamma0, epsB, top (quantile), tieTol
% Returns all diagnostics and the class R0 / R1 / D2 / I3 plus optimum relevance.
if nargin < 3 || isempty(opt), opt = struct(); end
d = @(f,v) getdef(opt,f,v);
delta = d('delta',0.01); gamma0 = d('gamma0',0.05); epsB = d('epsB',1e-9);
top = d('top',0.9); tieTol = d('tieTol',delta);   % designs within delta of the best are indistinguishable
P = P(:); B = B(:); A = P./B; n = numel(P);
r.n = n; r.validity = mean(P <= B*(1+1e-12));            % fraction with P <= B
r.CV_B = std(B)/mean(B); r.BRI = log(max(B)/min(B))/(2*delta);
rr = m05_rank_reversal(P, A, delta); r.Rrev = rr.Rrev; r.Rrev_delta = rr.Rrev_delta; r.tau = rr.tau;
r.rhoS = spearman_rho(P, A); r.pearson = corr1(P, A);
[~, o] = sort(P); dA = diff(A(o));
r.Gamma = sum(dA.^2)/(2*(n-1)*var(A) + eps);
r.NMI = nmi(P, A);
q = quantile1(P, top); it = P >= q;
if nnz(it) > 3
  rt = m05_rank_reversal(P(it), A(it), delta); r.Rrev_top = rt.Rrev_delta;
  Pt = P(it); At = A(it); [~, o2] = sort(Pt); dAt = diff(At(o2));
  r.Gamma_top = sum(dAt.^2)/(2*(numel(At)-1)*var(At) + eps);   % is A a function of P among the best designs?
else, r.Rrev_top = 0; r.Gamma_top = 0; end
% optimum-level regrets (optimistic over tie sets)
sP = P >= max(P)*(1 - tieTol); sA = A >= max(A)*(1 - tieTol);
r.dA_star = (max(A) - max(A(sP)))/max(A);     % attainment lost by maximising P
r.dP_star = (max(P) - max(P(sA)))/max(P);     % performance lost by maximising A
if nargin >= 4 && ~isempty(X)
  [~, iP] = max(P + 1e-12*A); [~, iA] = max(A + 1e-12*P);
  r.xP = X(iP,:); r.xA = X(iA,:);
end
% classification
if r.CV_B < epsB, r.class = 'R0';
elseif r.BRI <= 1 || r.Rrev_delta == 0, r.class = 'R1';
elseif r.Gamma < gamma0, r.class = 'D2';
else, r.class = 'I3'; end
r.shift = (r.dA_star > delta) || (r.dP_star > delta);
end
function v = getdef(s, f, v0)
if isfield(s, f), v = s.(f); else, v = v0; end
end
function q = quantile1(x, p)
x = sort(x(:)); q = x(max(1, min(numel(x), round(p*numel(x)))));
end
function c = corr1(a, b)
a = a - mean(a); b = b - mean(b); c = (a'*b)/sqrt((a'*a)*(b'*b) + eps);
end
function m = nmi(a, b)
n = numel(a); k = max(5, round(n^(1/3)));
ia = rankbin(a, k); ib = rankbin(b, k);
J = accumarray([ia ib], 1, [k k])/n; pa = sum(J,2); pb = sum(J,1);
nz = J > 0; PP = pa*pb; I = sum(J(nz).*log(J(nz)./PP(nz)));
Ha = -sum(pa(pa>0).*log(pa(pa>0))); Hb = -sum(pb(pb>0).*log(pb(pb>0)));
m = I/sqrt(Ha*Hb + eps);
end
function b = rankbin(x, k)
[~, o] = sort(x(:)); r = zeros(numel(x),1); r(o) = 1:numel(x);
b = min(k, 1 + floor((r-1)*k/numel(x)));
end
