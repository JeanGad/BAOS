function mask = nd_filter(F)
% ND_FILTER  Exact non-dominated mask (minimisation) for large sets.
% Points are visited in increasing order of their coordinate sum; a dominator
% always has a strictly smaller sum, so each point is checked only against the
% current front.  Exact duplicates are all kept.
n = size(F,1); [~, o] = sort(sum((F - min(F))./(max(F) - min(F) + eps), 2));
front = zeros(n,1); k = 0; mask = false(n,1);
for ii = 1:n
  i = o(ii); fi = F(i,:);
  if k > 0
    Fk = F(front(1:k),:);
    if any(all(Fk <= fi, 2) & any(Fk < fi, 2)), continue; end
  end
  k = k + 1; front(k) = i; mask(i) = true;
end
end
