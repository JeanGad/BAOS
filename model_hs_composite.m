function out = model_hs_composite(X, th)
% Two-phase isotropic composite (stiff phase 1, compliant phase 2).
% X = [f, theta]: f = stiff-phase volume fraction, theta = connectivity index
% of the stiff phase (reference shear modulus G0 = G2^(1-theta) * G1^theta).
% K* follows the Hashin-Shtrikman-type estimate with reference medium G0;
% theta = 0 and 1 reproduce the HS lower and upper bounds.
f = X(:,1); t = X(:,2);
hs = @(G0) 1./(f./(th.K1 + 4*G0/3) + (1-f)./(th.K2 + 4*G0/3)) - 4*G0/3;
G0 = th.G2.^(1-t).*th.G1.^t;
K = hs(G0); KU = hs(th.G1*ones(size(f))); KL = hs(th.G2*ones(size(f)));
rho = f*th.r1 + (1-f)*th.r2;                 % g/cm^3
out.P  = K./rho;                             % specific bulk modulus, GPa cm^3/g
out.BT = th.K1/th.r1*ones(size(f));          % theoretical bound (mediant inequality)
out.BA = KU./rho;                            % HS upper bound at the same f
out.BS = out.P./((K - KL)./(KU - KL));       % implied benchmark of bound-span attainment
out.K = K; out.KU = KU; out.KL = KL; out.rho = rho;
out.g = zeros(numel(f),0); out.cv = zeros(numel(f),1);
end
