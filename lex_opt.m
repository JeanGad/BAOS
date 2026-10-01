function [x, v] = lex_opt(ev, fp, fs, lb, ub, Xg, tol)
% LEX_OPT  Lexicographic optimum: maximise primary field fp, then maximise
% secondary field fs among designs whose primary value is within relative tol
% of the best. ev: X -> model output struct (with cv). Grid + Nelder-Mead.
o = ev(Xg); pen = 1e3*o.cv; fpg = -o.(fp) + pen; fpg(o.cv > 0) = inf;
f1 = @(x) pick(ev(x), fp, [], 0, 0);
[x1] = opt_single(f1, lb, ub, Xg, fpg); o1 = ev(x1); pstar = o1.(fp);
f2 = @(x) pick(ev(x), fs, fp, pstar, tol);
g2 = -o.(fs) + 1e3*max(pstar*(1 - tol) - o.(fp), 0)/abs(pstar) + pen; g2(o.cv > 0) = inf;
g2(o.(fp) < pstar*(1 - tol) - 1e-12) = inf;
if all(isinf(g2)), x = x1; else, x = opt_single(f2, lb, ub, Xg, g2); end
ox = ev(x); if ox.(fp) < pstar*(1 - tol) - 1e-9 || ox.cv > 0, x = x1; end
v = ev(x);
end
function f = pick(o, fs, fp, pstar, tol)
f = -o.(fs) + 1e3*o.cv;
if ~isempty(fp), f = f + 1e4*max(pstar*(1 - tol) - o.(fp), 0)/abs(pstar); end
end
