function prob = m01_define_problem(name)
% M01_DEFINE_PROBLEM  Replaceable {M, B, C, O} modules for each application.
% Returns: name, vars, lb, ub, th (parameters), eval (handle X->out),
% bench (list of benchmark field names), obj (objective catalogue with sense).
switch lower(name)
  case 'brayton'
    th = struct('T0',298.15,'p0',101.325,'cp',1.005,'k',1.4,'etac',0.85,'etat',0.88, ...
                'dph',0.03,'dpr0',0.01,'dpr1',0.03,'T4max',1000,'wmin',100,'EF',202);
    prob = struct('name','brayton','vars',{{'r_p','T_3 (K)','eps_r'}}, ...
                  'lb',[2 1100 0],'ub',[24 1700 0.95],'th',th);
    prob.eval = @(X) model_brayton(X, prob.th);
    prob.bench = {'BT','BA'}; prob.model = @model_brayton; prob.gscale = [th.T4max th.wmin]; prob.rfields = {'P','BT','BA','AT','AA','U','L'};
    prob.obj = struct('name',{{'P','AA','AT','U','L','C','E'}}, 'sense',[1 1 1 1 -1 -1 -1]);
  case 'hx'
    th = struct('mh',1.0,'cph',4.19,'cpc',4.18,'Thi',363.15,'Tci',293.15,'T0',293.15, ...
                'hh',3000,'hc0',3000,'mref',1.0,'Rw',1e-4,'dPc0',2e4,'dPh0',2e4, ...
                'Aref',20,'rho',1000,'etap',0.7);
    prob = struct('name','hx','vars',{{'m_c (kg/s)','A_hx (m^2)'}},'lb',[0.2 2],'ub',[3 40],'th',th);
    prob.eval = @(X) model_hx(X, prob.th);
    prob.bench = {'B'}; prob.model = @model_hx; prob.gscale = []; prob.rfields = {'P','B','A','U','L'};
    prob.obj = struct('name',{{'P','A','U','L','C'}}, 'sense',[1 1 1 -1 -1]);
  case 'solar'
    th = struct('Ac',2,'UL',5,'Fp',0.92,'ta',0.80,'G',800,'Ta',298.15,'cp',4180, ...
                'Tsun',5770,'Tomax',368.15);
    prob = struct('name','solar','vars',{{'m (kg/s)','T_i (K)'}},'lb',[0.005 298.15],'ub',[0.08 358.15],'th',th);
    prob.eval = @(X) model_solar(X, prob.th);
    prob.bench = {'B1','B2'}; prob.model = @model_solar; prob.gscale = th.Tomax; prob.rfields = {'P','B1','B2','A1','A2','U','L'};
    prob.obj = struct('name',{{'P','A1','A2','U','L'}}, 'sense',[1 1 1 1 -1]);
  otherwise
    error('unknown problem %s', name);
end
end
