function res = moead(fun, lb, ub, opt)
% MOEAD  MOEA/D with normalised Tchebycheff decomposition (Zhang & Li, 2007),
% SBX+PM operators, violation-first replacement, external non-dominated archive.
d = numel(lb); rng(opt.seed); M = opt.M; W = weights(M, opt.pop); N = size(W,1);
T = min(20, N); D = sqrt(max(sum(W.^2,2) + sum(W.^2,2)' - 2*(W*W'),0)); [~, B] = sort(D, 2); B = B(:,1:T);
X = lb + rand(N,d).*(ub - lb); [F, cv] = fun(X); z = min(F,[],1); nad = max(F,[],1);
res.trace = nan(opt.gen,1); res.nfe = N; AX = X; AF = F; AC = cv;
for g = 1:opt.gen
  Y = zeros(N,d); nb = zeros(N,T);
  for i = 1:N, k = B(i, randperm(T, 2)); Y(i,:) = variation(X(k(1),:), X(k(2),:), lb, ub, 1.0, 1/d, 20, 20); nb(i,:) = B(i, randperm(T)); end
  [FY, CY] = fun(Y); res.nfe = res.nfe + N;              % batch evaluation of offspring
  nad = max(AF,[],1);                                      % nadir from the non-dominated archive
  for i = 1:N
    fy = FY(i,:); cy = CY(i); z = min(z, fy); sc = max(nad - z, 1e-12);
    for j = nb(i,:)
      gy = max(W(j,:).*(fy - z)./sc); go = max(W(j,:).*(F(j,:) - z)./sc);
      if cy < cv(j) || (cy == cv(j) && gy <= go), X(j,:) = Y(i,:); F(j,:) = fy; cv(j) = cy; end
    end
  end
  AX = [AX; Y]; AF = [AF; FY]; AC = [AC; CY];
  [~, f1] = nd_sort(AF, AC); f1 = f1 & AC <= 0; AX = AX(f1,:); AF = AF(f1,:); AC = AC(f1);
  [~, u] = unique(round(AF*1e10)/1e10, 'rows'); AX = AX(u,:); AF = AF(u,:); AC = AC(u);
  while size(AF,1) > N, c = crowding(AF); [~, i] = min(c); AX(i,:) = []; AF(i,:) = []; AC(i) = []; end
  if isfield(opt,'trace') && ~isempty(opt.trace) && ~isempty(AF), res.trace(g) = opt.trace(AF); end
end
res.X = AX; res.F = AF;
end
function W = weights(M, N)
if M == 1, W = 1; return; end
H = 1; while nchoosek(H+M-1, M-1) < N, H = H + 1; end
C = nchoosek(1:H+M-1, M-1); W = zeros(size(C,1), M); W(:,1) = C(:,1) - 1;
for m = 2:M-1, W(:,m) = C(:,m) - C(:,m-1) - 1; end
W(:,M) = H + M - 1 - C(:,M-1); W = W/H; W(W == 0) = 1e-6;
end
