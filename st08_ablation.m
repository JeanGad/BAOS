function Ab = st08_ablation(cfg)
% ST08_ABLATION  Objective-set ablation (Brayton, exact grid Pareto sets):
% A: P | B: P+A | C: P+A+L | D: P+A+L+U | E: P+A+L+U+C, for A in {A_T, A_A},
% each compared with the same set WITHOUT the benchmark objective.
p = m01_define_problem('brayton'); Xg = m07_doe(31^3, p.lb, p.ub, 'grid', 1); o = p.eval(Xg);
f = o.cv <= 0; Xg = Xg(f,:); o = p.eval(Xg);
base = {{}, {'L'}, {'L','U'}, {'L','U','C'}}; sgn = struct('P',-1,'AT',-1,'AA',-1,'U',-1,'L',1,'C',1);
rows = [];
for b = {'AT','AA'}
  for k = 1:numel(base)
    wo = [{'P'} base{k}]; wi = [{'P', b{1}} base{k}];
    FW = mk(o, wo, sgn); FB = mk(o, wi, sgn); mW = nd_filter(FW); mB = nd_filter(FB);
    tw = m14_mcdm(-FW(mW,:), true(1,numel(wo)), ones(1,numel(wo))/numel(wo));
    tb = m14_mcdm(-FB(mB,:), true(1,numel(wi)), ones(1,numel(wi))/numel(wi));
    XW = Xg(mW,:); XB = Xg(mB,:); xw = XW(tw.best,:); xb = XB(tb.best,:); ow = p.eval(xw); ob = p.eval(xb);
    rows = [rows; strcmp(b{1},'AA')+1, k+1, nnz(mW), nnz(mB), nnz(mW & mB)/nnz(mW | mB), nnz(mB & ~mW)/nnz(mB), ...
        norm((xb - xw)./(p.ub - p.lb))/sqrt(3), ow.P, ob.P, ow.(b{1}), ob.(b{1}), xw, xb];
  end
end
write_csv(fullfile(cfg.out, 'st08_ablation.csv'), {'bench(1=AT,2=AA)','set(2=P+A,3=+L,4=+U,5=+C)','nPareto_without','nPareto_with','O_P','phi_B', ...
   'D_compromise','P_comp_without','P_comp_with','A_comp_without','A_comp_with','xw_rp','xw_T3','xw_er','xb_rp','xb_T3','xb_er'}, rows);
Ab.rows = rows;
end
function F = mk(o, names, sgn)
F = zeros(numel(o.P), numel(names)); for j = 1:numel(names), F(:,j) = sgn.(names{j})*o.(names{j}); end
end
