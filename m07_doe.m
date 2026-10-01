function X = m07_doe(n, lb, ub, method, seed)
% M07_DOE  Space-filling designs: 'lhs' (maximin over 25 candidates),
% 'halton' (low-discrepancy, used in place of Sobol sequences), 'grid', 'random'.
d = numel(lb); rng(seed);
switch method
  case 'lhs'
    best = -inf;
    for c = 1:25
      U = zeros(n,d);
      for j = 1:d, U(:,j) = (randperm(n)' - rand(n,1))/n; end
      D = pdist2s(U); D(1:n+1:end) = inf; md = min(D(:));
      if md > best, best = md; Ub = U; end
    end
    U = Ub;
  case 'halton'
    pr = [2 3 5 7 11 13 17 19]; U = zeros(n,d); skip = 20;
    for j = 1:d, U(:,j) = vdc((1:n)' + skip, pr(j)); end
  case 'grid'
    k = round(n^(1/d)); g = linspace(0,1,k); C = cell(1,d); [C{:}] = ndgrid(g);
    U = zeros(numel(C{1}), d); for j = 1:d, U(:,j) = C{j}(:); end
  otherwise
    U = rand(n,d);
end
X = lb + U.*(ub - lb);
end
function v = vdc(i, b)
v = zeros(size(i)); f = 1/b; i = i(:);
while any(i > 0), v = v + f*mod(i,b); i = floor(i/b); f = f/b; end
end
function D = pdist2s(U)
s = sum(U.^2,2); D = sqrt(max(s + s' - 2*(U*U'), 0));
end
