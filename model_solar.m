function out = model_solar(X, th)
% MODEL_SOLAR  Flat-plate collector, Hottel-Whillier-Bliss.
% X = [m (kg/s), T_i (K)]. P = collector efficiency.
% B1 = tau*alpha (optical, zero-loss limit: constant  -> expected redundancy)
% B2 = F_R*tau*alpha (zero-loss limit at the SAME flow rate: design-dependent)
m = X(:,1); Ti = X(:,2); cp = th.cp;
FR = m*cp./(th.Ac*th.UL).*(1 - exp(-th.Ac*th.UL*th.Fp./(m*cp)));
Qu = th.Ac.*FR.*(th.G*th.ta - th.UL.*(Ti - th.Ta));
eta = Qu./(th.Ac*th.G); To = Ti + Qu./(m*cp); T0 = th.Ta;
exU = m*cp.*((To - Ti) - T0.*log(To./Ti));
psiP = 1 - 4/3*(T0/th.Tsun) + 1/3*(T0/th.Tsun)^4;              % Petela factor
out.P = eta; out.B1 = th.ta*ones(size(m)); out.B2 = FR*th.ta;
out.A1 = out.P./out.B1; out.A2 = out.P./out.B2;
out.U = exU/1e3; out.L = (th.Ac*th.G*psiP - exU)/1e3; out.To = To; out.FR = FR;
out.g = To - th.Tomax; out.cv = max(out.g,0)/th.Tomax;
end
