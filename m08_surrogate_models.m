function varargout = m08_surrogate_models(action, varargin)
% M08_SURROGATE_MODELS  Toolbox-free surrogate library.
%   mdl = m08_surrogate_models('fit', type, X, y, lb, ub, seed)
%   [mu, sd] = m08_surrogate_models('predict', mdl, X)
% types: 'RSM' (full quadratic), 'GPR' (ARD squared-exponential, ML-II),
%        'ANN' (1 hidden tanh layer, Adam, 3-member ensemble),
%        'RF' (bagged CART, 200 trees), 'GBM' (gradient-boosted depth-3 trees).
% SVR is NOT implemented (would need a QP solver); see limitations.
switch action
  case 'fit'
    [type, X, y, lb, ub, seed] = varargin{:}; rng(seed);
    Z = (X - lb)./(ub - lb); my = mean(y); sy = std(y) + eps; t = (y - my)/sy;
    mdl = struct('type', type, 'lb', lb, 'ub', ub, 'my', my, 'sy', sy);
    switch type
      case 'RSM', Phi = quad(Z); mdl.w = (Phi'*Phi + 1e-8*eye(size(Phi,2)))\(Phi'*t);
      case 'GPR', mdl = gpr_fit(mdl, Z, t);
      case 'ANN', mdl.nets = cell(1,3); for e = 1:3, mdl.nets{e} = ann_fit(Z, t, 12, 2500, seed + e); end
      case 'RF'
        mdl.trees = cell(1,200); n = size(Z,1);
        for b = 1:200, i = randi(n, n, 1); mdl.trees{b} = tree_fit(Z(i,:), t(i), 3, 30, max(1, size(Z,2)-1)); end
      case 'GBM'
        mdl.f0 = mean(t); r = t - mdl.f0; mdl.lr = 0.05; mdl.trees = cell(1,400);
        for b = 1:400
          i = randperm(size(Z,1), round(0.8*size(Z,1)));
          tr = tree_fit(Z(i,:), r(i), 3, 3, size(Z,2)); mdl.trees{b} = tr; r = r - mdl.lr*tree_pred(tr, Z);
        end
    end
    varargout{1} = mdl;
  case 'predict'
    [mdl, X] = varargin{:}; Z = (X - mdl.lb)./(mdl.ub - mdl.lb); sd = nan(size(X,1),1);
    switch mdl.type
      case 'RSM', t = quad(Z)*mdl.w;
      case 'GPR', [t, v] = gpr_pred(mdl, Z); sd = sqrt(max(v,0))*mdl.sy;
      case 'ANN', T = zeros(size(Z,1),3); for e = 1:3, T(:,e) = ann_pred(mdl.nets{e}, Z); end; t = mean(T,2); sd = std(T,0,2)*mdl.sy;
      case 'RF', T = zeros(size(Z,1),numel(mdl.trees)); for b = 1:numel(mdl.trees), T(:,b) = tree_pred(mdl.trees{b}, Z); end; t = mean(T,2); sd = std(T,0,2)*mdl.sy;
      case 'GBM', t = mdl.f0*ones(size(Z,1),1); for b = 1:numel(mdl.trees), t = t + mdl.lr*tree_pred(mdl.trees{b}, Z); end
    end
    varargout{1} = mdl.my + mdl.sy*t; varargout{2} = sd;
end
end
% ---------------- RSM ----------------
function Phi = quad(Z)
[n, d] = size(Z); Phi = [ones(n,1) Z Z.^2];
for i = 1:d-1, for j = i+1:d, Phi = [Phi Z(:,i).*Z(:,j)]; end, end
end
% ---------------- GPR ----------------
function mdl = gpr_fit(mdl, Z, t)
d = size(Z,2); best = inf; f = @(h) gp_nll(h, Z, t);
starts = [log(0.3)*ones(1,d) 0 log(1e-3); log(1)*ones(1,d) 0 log(1e-2); log(0.1)*ones(1,d) 0 log(1e-4)];
o = optimset('MaxFunEvals', 1500, 'MaxIter', 1500, 'TolX', 1e-6, 'TolFun', 1e-8, 'Display', 'off');
for s = 1:size(starts,1)
  [h, v] = fminsearch(f, starts(s,:), o); if v < best, best = v; hb = h; end
end
mdl.h = hb; mdl.Z = Z; [~, mdl.L, mdl.alpha] = gp_nll(hb, Z, t); mdl.nll = best;
end
function [v, L, alpha] = gp_nll(h, Z, t)
K = kern(Z, Z, h); n = size(Z,1); sn2 = exp(2*h(end)) + 1e-8;
[L, p] = chol(K + sn2*eye(n), 'lower'); if p > 0, v = 1e10; alpha = []; return; end
alpha = L'\(L\t); v = 0.5*t'*alpha + sum(log(diag(L))) + 0.5*n*log(2*pi);
if any(abs(h(1:end-2)) > 7), v = v + 1e3; end            % keep length-scales in a sane range
end
function K = kern(A, B, h)
d = size(A,2); l = exp(h(1:d)); sf2 = exp(2*h(d+1));
A = A./l; B = B./l; D = sum(A.^2,2) + sum(B.^2,2)' - 2*A*B'; K = sf2*exp(-0.5*max(D,0));
end
function [m, v] = gpr_pred(mdl, Z)
Ks = kern(Z, mdl.Z, mdl.h); m = Ks*mdl.alpha; W = mdl.L\Ks';
v = exp(2*mdl.h(numel(mdl.h)-1)) - sum(W.^2,1)';
end
% ---------------- ANN ----------------
function net = ann_fit(Z, t, H, iters, seed)
rng(seed); [n, d] = size(Z); nv = max(3, round(0.15*n)); p = randperm(n); iv = p(1:nv); it = p(nv+1:end);
W1 = randn(d,H)*sqrt(1/d); b1 = zeros(1,H); W2 = randn(H,1)*sqrt(1/H); b2 = 0;
th = {W1,b1,W2,b2}; m = cellfun(@(x) 0*x, th, 'UniformOutput', false); v = m;
lr = 0.01; b1a = 0.9; b2a = 0.999; lam = 1e-4; best = inf; bestth = th;
for k = 1:iters
  [g, ~] = ann_grad(th, Z(it,:), t(it), lam);
  for q = 1:4
    m{q} = b1a*m{q} + (1-b1a)*g{q}; v{q} = b2a*v{q} + (1-b2a)*g{q}.^2;
    th{q} = th{q} - lr*(m{q}/(1-b1a^k))./(sqrt(v{q}/(1-b2a^k)) + 1e-8);
  end
  if mod(k,25) == 0
    ev = mean((ann_fwd(th, Z(iv,:)) - t(iv)).^2); if ev < best, best = ev; bestth = th; end
  end
end
net = bestth;
end
function y = ann_fwd(th, Z)
y = tanh(Z*th{1} + th{2})*th{3} + th{4};
end
function [g, L] = ann_grad(th, Z, t, lam)
n = size(Z,1); A1 = tanh(Z*th{1} + th{2}); y = A1*th{3} + th{4}; e = (y - t)/n; L = sum((y-t).^2)/n;
g{4} = 2*sum(e); g{3} = 2*A1'*e + 2*lam*th{3}; dA = (2*e*th{3}').*(1 - A1.^2);
g{2} = sum(dA,1); g{1} = Z'*dA + 2*lam*th{1};
end
function y = ann_pred(net, Z)
y = ann_fwd(net, Z);
end
% ---------------- CART (used by RF and GBM) ----------------
function T = tree_fit(Z, t, minleaf, maxdepth, mtry)
[n, d] = size(Z); cap = 2*n; T.f = zeros(cap,1); T.s = zeros(cap,1); T.l = zeros(cap,1); T.r = zeros(cap,1); T.v = zeros(cap,1);
stack = {1:n}; depth = 0; nodes = 1; q = {1}; dep = 0; id = 1;
st_idx = {1:n}; st_node = 1; st_dep = 0; top = 1;
while top > 0
  idx = st_idx{top}; node = st_node(top); dd = st_dep(top); top = top - 1;
  y = t(idx); T.v(node) = mean(y);
  if numel(idx) < 2*minleaf || dd >= maxdepth || var(y) < 1e-14, continue; end
  best = inf; feats = randperm(d, mtry);
  for f = feats
    [zs, o] = sort(Z(idx,f)); ys = y(o); cs = cumsum(ys); cs2 = cumsum(ys.^2); m = numel(ys);
    k = (minleaf:m-minleaf)'; if isempty(k), continue; end
    sl = cs2(k) - cs(k).^2./k; sr = (cs2(m) - cs2(k)) - (cs(m) - cs(k)).^2./(m - k);
    ok = zs(k) < zs(k+1); sse = sl + sr; sse(~ok) = inf; [v, j] = min(sse);
    if v < best, best = v; bf = f; bs = 0.5*(zs(k(j)) + zs(k(j)+1)); end
  end
  if isinf(best), continue; end
  L = idx(Z(idx,bf) <= bs); R = idx(Z(idx,bf) > bs);
  T.f(node) = bf; T.s(node) = bs; T.l(node) = id + 1; T.r(node) = id + 2;
  top = top + 1; st_idx{top} = L; st_node(top) = id + 1; st_dep(top) = dd + 1;
  top = top + 1; st_idx{top} = R; st_node(top) = id + 2; st_dep(top) = dd + 1; id = id + 2;
end
T.f = T.f(1:id); T.s = T.s(1:id); T.l = T.l(1:id); T.r = T.r(1:id); T.v = T.v(1:id);
end
function y = tree_pred(T, Z)
n = size(Z,1); node = ones(n,1); act = T.f(node) > 0;
while any(act)
  a = find(act); f = T.f(node(a)); zz = Z(sub2ind(size(Z), a, f));
  goL = zz <= T.s(node(a)); node(a(goL)) = T.l(node(a(goL))); node(a(~goL)) = T.r(node(a(~goL)));
  act = T.f(node) > 0;
end
y = T.v(node);
end
