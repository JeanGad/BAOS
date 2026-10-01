function T = m18_statistics(samples, names)
% M18_STATISTICS  Pairwise two-sided Mann-Whitney U tests (normal approximation
% with tie correction), Holm-Bonferroni adjustment and Vargha-Delaney A12.
% samples: runs x k matrix (one column per algorithm).
k = size(samples,2); T = struct('a',{},'b',{},'p',{},'padj',{},'A12',{});
ps = [];
for i = 1:k-1, for j = i+1:k
  x = samples(:,i); y = samples(:,j); nx = numel(x); ny = numel(y);
  r = tiedrank1([x; y]); U = sum(r(1:nx)) - nx*(nx+1)/2; mu = nx*ny/2;
  [~,~,g] = unique([x; y]); t = accumarray(g,1); N = nx + ny;
  sig = sqrt(nx*ny/12*((N+1) - sum(t.^3 - t)/(N*(N-1))));
  z = (U - mu - 0.5*sign(U - mu))/max(sig, eps); p = erfc(abs(z)/sqrt(2));
  T(end+1) = struct('a',names{i},'b',names{j},'p',p,'padj',NaN,'A12',U/(nx*ny)); ps(end+1) = p;
end, end
[ps2, o] = sort(ps); m = numel(ps); adj = zeros(1,m); run = 0;
for q = 1:m, run = max(run, min(1, (m - q + 1)*ps2(q))); adj(o(q)) = run; end
for q = 1:m, T(q).padj = adj(q); end
end
