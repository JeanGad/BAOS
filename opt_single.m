function [xb, fb] = opt_single(f, lb, ub, Xg, fg)
% OPT_SINGLE  Deterministic single-objective optimum: best of a dense grid
% (f already evaluated as fg, minimisation, infeasible = +inf) polished by
% Nelder-Mead in normalised coordinates with bound clipping.
[~, i] = min(fg); z0 = (Xg(i,:) - lb)./(ub - lb); fr = ub > lb;
g = @(z) f(lb + min(max(expand(z, z0, fr),0),1).*(ub - lb));
o = optimset('TolX',1e-10,'TolFun',1e-12,'MaxFunEvals',4000,'MaxIter',4000,'Display','off');
zb = fminsearch(g, z0(fr), o); zz = min(max(expand(zb, z0, fr),0),1);
xb = lb + zz.*(ub - lb); fb = f(xb);
if fb > fg(i), xb = Xg(i,:); fb = fg(i); end
end
function z = expand(zf, z0, fr)
z = z0; z(fr) = zf;
end
