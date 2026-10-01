function C = variation(P1, P2, lb, ub, pc, pm, etac, etam)
% VARIATION  SBX crossover + polynomial mutation (vectorised, bounded).
[n, d] = size(P1); u = rand(n,d);
beta = (2*u).^(1/(etac+1)); hi = u > 0.5; beta(hi) = (1./(2*(1-u(hi)))).^(1/(etac+1));
beta(rand(n,d) > 0.5) = 1; beta(repmat(rand(n,1) > pc, 1, d)) = 1;
C = 0.5*((1+beta).*P1 + (1-beta).*P2);
mut = rand(n,d) < pm; u = rand(n,d); R = repmat(ub - lb, n, 1);
Rs = max(ub - lb, eps); dl = (C - lb)./Rs; du = (ub - C)./Rs; dq = zeros(n,d); lo = u < 0.5;
dq(lo)  = (2*u(lo) + (1-2*u(lo)).*(1-dl(lo)).^(etam+1)).^(1/(etam+1)) - 1;
dq(~lo) = 1 - (2*(1-u(~lo)) + 2*(u(~lo)-0.5).*(1-du(~lo)).^(etam+1)).^(1/(etam+1));
C(mut) = C(mut) + dq(mut).*R(mut); C = min(max(C, lb), ub);
end
