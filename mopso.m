function res = mopso(fun, lb, ub, opt)
% MOPSO  Multi-objective PSO with crowding-distance archive and leader
% selection (Raquel & Naval, 2005) plus decaying polynomial-mutation turbulence.
d = numel(lb); N = opt.pop; rng(opt.seed); Na = N;
X = lb + rand(N,d).*(ub - lb); V = zeros(N,d); [F, cv] = fun(X);
PX = X; PF = F; PC = cv; res.nfe = N; res.trace = nan(opt.gen,1);
[~, f1] = nd_sort(F, cv); k = f1 & cv <= 0;
if any(k), AX = X(k,:); AF = F(k,:); else, [~, i] = min(cv); AX = X(i,:); AF = F(i,:); end
for g = 1:opt.gen
  w = 0.9 - 0.5*g/opt.gen; cd = crowding(AF); fin = ~isinf(cd);
  if any(fin), cd(~fin) = 2*max(cd(fin)); else, cd(:) = 1; end
  cs = cumsum(cd/sum(cd)); L = zeros(N,1);
  for i = 1:N, L(i) = find(cs >= rand*cs(end), 1); end
  V = w*V + 1.5*rand(N,d).*(PX - X) + 1.5*rand(N,d).*(AX(L,:) - X);
  vmax = 0.5*(ub - lb); V = max(min(V, vmax), -vmax); X = X + V;
  out = X < lb | X > ub; V(out) = -V(out); X = min(max(X, lb), ub);
  mi = rand(N,1) < 0.5*(1 - (g-1)/opt.gen)^1.5;
  if any(mi), X(mi,:) = variation(X(mi,:), X(mi,:), lb, ub, 0, 1/d, 20, 20); end
  [F, cv] = fun(X); res.nfe = res.nfe + N;
  bt = (cv < PC) | (cv == PC & all(F <= PF, 2) & any(F < PF, 2));
  nd = cv == PC & ~(all(PF <= F, 2) & any(PF < F, 2)) & ~bt & rand(N,1) < 0.5;
  up = bt | nd; PX(up,:) = X(up,:); PF(up,:) = F(up,:); PC(up) = cv(up);
  fe = cv <= 0; AX = [AX; X(fe,:)]; AF = [AF; F(fe,:)];
  [~, f1] = nd_sort(AF); AX = AX(f1,:); AF = AF(f1,:);
  [~, u] = unique(round(AF*1e10)/1e10, 'rows'); AX = AX(u,:); AF = AF(u,:);
  while size(AF,1) > Na, c = crowding(AF); [~, i] = min(c); AX(i,:) = []; AF(i,:) = []; end
  if isfield(opt,'trace') && ~isempty(opt.trace), res.trace(g) = opt.trace(AF); end
end
res.X = AX; res.F = AF;
end
