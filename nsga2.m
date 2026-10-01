function res = nsga2(fun, lb, ub, opt)
% NSGA2  Elitist non-dominated sorting GA (Deb et al., 2002), toolbox-free.
% fun: X -> [F (minimise), cv]. opt: pop, gen, seed, pc, pm, etac, etam, trace (handle F->scalar)
N = opt.pop; d = numel(lb); rng(opt.seed);
pc = getf(opt,'pc',0.9); pm = getf(opt,'pm',1/d); ec = getf(opt,'etac',15); em = getf(opt,'etam',20);
X = lb + rand(N,d).*(ub - lb); [F, cv] = fun(X); res.trace = nan(opt.gen,1); res.nfe = N;
rk = nd_sort(F, cv); cd = allcrowd(F, rk);
for g = 1:opt.gen
  i1 = tourn(rk, cd, N); i2 = tourn(rk, cd, N);
  Q = variation(X(i1,:), X(i2,:), lb, ub, pc, pm, ec, em);
  [FQ, cvQ] = fun(Q); res.nfe = res.nfe + N;
  XX = [X; Q]; FF = [F; FQ]; CC = [cv; cvQ];
  rk2 = nd_sort(FF, CC); cd2 = allcrowd(FF, rk2);
  [~, o] = sortrows([rk2, -cd2]); o = o(1:N);
  X = XX(o,:); F = FF(o,:); cv = CC(o); rk = rk2(o); cd = cd2(o);
  if isfield(opt,'trace') && ~isempty(opt.trace)
    m = rk == 1 & cv <= 0; if any(m), res.trace(g) = opt.trace(F(m,:)); end
  end
end
m = rk == 1 & cv <= 0; res.X = X(m,:); res.F = F(m,:);
[~, u] = unique(round(res.F*1e10)/1e10, 'rows'); res.X = res.X(u,:); res.F = res.F(u,:);
end
function cd = allcrowd(F, rk)
cd = zeros(size(F,1),1); for r = unique(rk)', m = rk == r; cd(m) = crowding(F(m,:)); end
end
function w = tourn(rk, cd, N)
a = randi(numel(rk), N, 1); b = randi(numel(rk), N, 1);
better = rk(a) < rk(b) | (rk(a) == rk(b) & cd(a) > cd(b)); w = b; w(better) = a(better);
end
function v = getf(s, f, v0)
if isfield(s,f), v = s.(f); else, v = v0; end
end
