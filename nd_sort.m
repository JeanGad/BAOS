function [rank, front1] = nd_sort(F, cv)
% ND_SORT  Non-dominated sorting with Deb's constraint-domination.
% F: n x M (minimise), cv: n x 1 total violation (0 = feasible).
n = size(F,1); if nargin < 2, cv = zeros(n,1); end
D = false(n);                                     % D(i,j): i dominates j
for i = 1:n
  le = all(F(i,:) <= F, 2); lt = any(F(i,:) < F, 2); D(i,:) = (le & lt)';
end
feas = cv <= 0; D(feas, ~feas) = true; D(~feas, feas) = false;
inf2 = find(~feas);
if ~isempty(inf2), C = cv(inf2); D(inf2, inf2) = C < C'; end
rank = zeros(n,1); cnt = sum(D,1)'; r = 1; cur = find(cnt == 0);
while ~isempty(cur)
  rank(cur) = r; cnt = cnt - sum(D(cur,:),1)'; cnt(rank > 0) = -1;
  cur = find(cnt == 0); r = r + 1;
end
front1 = rank == 1;
end
