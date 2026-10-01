function v = hv(F, r)
% HV  Hypervolume (minimisation) w.r.t. reference r. Exact for M<=2 and for M=3 with <=400 points,
% Monte Carlo (5e4 points, fixed seed, chunked) for M>3.
F = F(all(F < r, 2), :); M = numel(r);
if isempty(F), v = 0; return; end
if M == 1, v = r - min(F); return; end
if M == 2, v = hv2(F, r); return; end
if M == 3 && size(F,1) <= 400
  [~, o] = sort(F(:,3)); F = F(o,:); z = [F(:,3); r(3)]; v = 0;
  for i = 1:size(F,1), v = v + hv2(F(1:i,1:2), r(1:2))*(z(i+1) - z(i)); end
  return
end
s = rng; rng(12345); K = 5e4; lo = min(F,[],1); S = lo + rand(K,M).*(r - lo); dom = false(K,1);
Ft = permute(F, [3 1 2]);                                  % 1 x nF x M
for c = 1:1000:K
  ii = c:min(K, c+999); Sc = permute(S(ii,:), [1 3 2]); % nc x 1 x M
  dom(ii) = any(all(bsxfun(@ge, Sc, Ft), 3), 2);
end
v = mean(dom)*prod(r - lo); rng(s);
end
function v = hv2(F, r)
[~, o] = sort(F(:,1)); F = F(o,:); v = 0; best = r(2);
for i = 1:size(F,1)
  if F(i,2) < best, v = v + (r(1) - F(i,1))*(best - F(i,2)); best = F(i,2); end
end
end
