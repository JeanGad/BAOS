function out = model_hx(X, th)
% MODEL_HX  Counterflow liquid-liquid heat exchanger (epsilon-NTU), cold-side
% flow-dependent film coefficient, friction losses scaling with length.
% X = [m_c (kg/s), A_hx (m2)]. P = duty Q (kW); B = Q_max = C_min*dT_max
% (infinite-area limit at the SAME flows, i.e. conditionally attainable).
mc = X(:,1); Ahx = X(:,2);
Ch = th.mh*th.cph*1e3*ones(size(mc)); Cc = mc*th.cpc*1e3;       % W/K
hc = th.hc0*(mc/th.mref).^0.8;  U = 1./(1/th.hh + 1./hc + th.Rw);
UA = U.*Ahx;  Cmin = min(Ch,Cc); Cmax = max(Ch,Cc); Cr = Cmin./Cmax; NTU = UA./Cmin;
e = zeros(size(mc)); one = abs(1-Cr) < 1e-9;
e(~one) = (1 - exp(-NTU(~one).*(1-Cr(~one))))./(1 - Cr(~one).*exp(-NTU(~one).*(1-Cr(~one))));
e(one)  = NTU(one)./(1+NTU(one));
Qmax = Cmin.*(th.Thi - th.Tci); Q = e.*Qmax;
Tho = th.Thi - Q./Ch; Tco = th.Tci + Q./Cc; T0 = th.T0;
SgenHT = Ch.*log(Tho./th.Thi) + Cc.*log(Tco./th.Tci);          % W/K
dPc = th.dPc0*(mc/th.mref).^1.8.*(Ahx/th.Aref); dPh = th.dPh0*(Ahx/th.Aref);
Wp  = (mc.*dPc + th.mh.*dPh)./(th.rho*th.etap);                 % W electric
exGain = Cc.*((Tco - th.Tci) - T0.*log(Tco./th.Tci));           % exergy to cold stream
exRel  = Ch.*((th.Thi - Tho) - T0.*log(th.Thi./Tho));           % exergy released by hot stream
out.P = Q/1e3; out.B = Qmax/1e3; out.A = out.P./out.B;
out.U = exGain/1e3; out.L = (T0.*SgenHT + Wp)/1e3; out.C = Ahx; out.Wp = Wp/1e3;
out.Sgen = SgenHT; out.closure = (exRel - exGain - T0.*SgenHT)/1e3;
out.parts = [T0.*SgenHT, Wp]/1e3; out.partnames = {'heat-transfer exergy destruction','pumping'};
out.g = zeros(numel(mc),0); out.cv = zeros(size(mc));
end
