function s = m16_optimum_shift(Xg, FC, FB, lb, ub, xP, xA, PA)
% M16_OPTIMUM_SHIFT  Optimum-shift metrics on a common evaluated design set Xg.
% FC: conventional objective matrix (minimise), FB: benchmark-aware (minimise).
% xP, xA: single-objective optima of P and A. PA = [P A] at xP and xA rows.
mC = nd_filter(FC); mB = nd_filter(FB);
s.nC = nnz(mC); s.nB = nnz(mB);
s.O_P = nnz(mC & mB)/nnz(mC | mB);                      % Jaccard overlap of Pareto sets
s.phi_B = nnz(mB & ~mC)/nnz(mB);                        % benchmark-induced fraction
s.OSI = norm((xA - xP)./(ub - lb))/sqrt(numel(lb));     % normalised to [0,1]
s.D_x = norm((xA - xP)./(ub - lb));
s.dP = (PA(1,1) - PA(2,1))/PA(1,1); s.dA = (PA(2,2) - PA(1,2))/PA(2,2);
% hypervolume of both sets in the benchmark-aware space (common normalisation)
Z = FB([find(mC); find(mB)],:); id = min(Z,[],1); nd = max(Z,[],1);
s.HV_C = hv((FB(mC,:) - id)./(nd - id + eps), 1.1*ones(1,size(FB,2)));
s.HV_B = hv((FB(mB,:) - id)./(nd - id + eps), 1.1*ones(1,size(FB,2)));
s.dHV = s.HV_B - s.HV_C; s.maskC = mC; s.maskB = mB;
end
