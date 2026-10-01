function out = model_brayton(X, th)
% MODEL_BRAYTON  Regenerative Brayton cycle, cold-air-standard, with component
% isentropic efficiencies and pressure losses. Exact exergy bookkeeping
% (Gouy-Stodola) so that B_A - eta is fully attributable.
% X = [rp, T3 (K), eps_r]  (rows = designs).  th = parameter struct.
% Outputs are column vectors (kJ/kg for specific quantities).
rp = X(:,1); T3 = X(:,2); er = X(:,3);
T0 = th.T0; cp = th.cp; k = th.k; R = cp*(k-1)/k; a = (k-1)/k;
T1 = T0; p0 = th.p0;
T2s = T1.*rp.^a;             T2 = T1 + (T2s - T1)./th.etac;
dpr = th.dpr0 + th.dpr1.*er.^2;          % regenerator pressure-loss fraction (each side)
p2 = p0.*rp;  px = p2.*(1-dpr);  p3 = px.*(1-th.dph);
p4 = p0./(1-dpr);                          % turbine back-pressure (hot-side regenerator loss)
T4s = T3.*(p4./p3).^a;       T4 = T3 - th.etat.*(T3 - T4s);
dT = max(T4 - T2, 0);        Tx = T2 + er.*dT;      T5 = T4 - (Tx - T2);
qin = cp.*(T3 - Tx);  wc = cp.*(T2 - T1);  wt = cp.*(T3 - T4);  wnet = wt - wc;
eta = wnet./qin;
s = @(T,p) cp.*log(T./T0) - R.*log(p./p0);   % specific entropy rel. to dead state
s1 = 0; s2 = s(T2,p2); sx = s(Tx,px); s3 = s(T3,p3); s4 = s(T4,p4); s5 = s(T5,p0);
edc = T0.*(s2 - s1);  edt = T0.*(s4 - s3);  edr = T0.*((sx - s2) + (s5 - s4));
edh = T0.*(s3 - sx - cp.*log(T3./Tx));        % heater: pressure-loss part only
exh = cp.*(T5 - T0) - T0.*cp.*log(T5./T0);   % exergy rejected with exhaust
exq = qin - T0.*cp.*log(T3./Tx);              % exergy of heat actually added
out.P   = eta;
out.BT  = 1 - T0./T3;                          % theoretical Carnot bound (T3 is a decision variable)
out.BA  = exq./qin;                            % conditionally attainable limit (Carnot factor of actual heat input)
out.AT  = out.P./out.BT;  out.AA = out.P./out.BA;
out.U   = wnet;  out.L = edc + edt + edr + edh;   out.Lex = out.L + exh;
out.C   = qin;   out.E = th.EF./max(eta,1e-6);     % g CO2 / kWh_e (strictly decreasing in eta)
out.parts = [out.BT - out.BA, edc./qin, edt./qin, edr./qin, edh./qin, exh./qin]; % efficiency points
out.partnames = {'context (B_T-B_A)','compressor','turbine','regenerator','heater dP','exhaust'};
out.closure = exq - wnet - (edc + edt + edr + edh + exh);
out.T4 = T4; out.T2 = T2; out.Tx = Tx; out.T5 = T5; out.Sgen = (edc+edt+edr+edh)./T0;
out.g = [T4 - th.T4max, th.wmin - wnet];      % g <= 0 feasible
out.cv = sum(max(out.g,0)./[th.T4max, th.wmin], 2);
end
