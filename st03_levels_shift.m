function L = st03_levels_shift(cfg, E)
% ST03_LEVELS_SHIFT  Objective selection (m06), three optimisation levels,
% optimum-shift metrics (m16), MCDM compromise (m14), gap attribution and the
% engineering validation table, for each engineering model.
spec = { 'brayton', {'P','AT','AA','U','L','C','E'}, [1 1 1 1 -1 -1 -1], {'AT','AA'}, 31, [1 2];
         'hx',      {'P','A','U','L','C'},            [1 1 1 -1 -1],      {'A'},       121, 5;
         'solar',   {'P','A1','A2','U','L'},          [1 1 1 1 -1],       {'A1','A2'}, 121, [6 7] };
L.sel = {}; L.shift = []; L.shiftlab = {}; L.table = []; L.tablab = {};
for m = 1:size(spec,1)
  p = m01_define_problem(spec{m,1}); names = spec{m,2}; sense = spec{m,3}; printf('  st03: %s\n', p.name); fflush(stdout);
  rng(cfg.seed); X = p.lb + rand(4*cfg.nBOIT, numel(p.lb)).*(p.ub - p.lb); o = p.eval(X); X = X(o.cv <= 0,:); X = X(1:cfg.nBOIT,:); o = p.eval(X);
  Y = zeros(size(X,1), numel(names)); for j = 1:numel(names), Y(:,j) = o.(names{j}); end
  sel = m06_objective_selection(Y, names, sense, struct('delta', cfg.delta)); L.sel{m} = sel;
  write_csv(fullfile(cfg.out, sprintf('st03_selection_%s_Rrevdelta.csv', p.name)), names, sel.Rrev_delta);
  write_csv(fullfile(cfg.out, sprintf('st03_selection_%s_rhoS.csv', p.name)), names, sel.rhoS);
  kept = names(sel.keep); ks = sense(sel.keep);
  conv = ~ismember(kept, spec{m,4});                        % conventional objective set
  Xg = m07_doe(spec{m,5}^numel(p.lb), p.lb, p.ub, 'grid', 1); og = p.eval(Xg); fe = og.cv <= 0; Xg = Xg(fe,:); og = p.eval(Xg);
  FC = zeros(size(Xg,1), nnz(conv)); cc = find(conv); for j = 1:numel(cc), FC(:,j) = -ks(cc(j))*og.(kept{cc(j)}); end
  for bb = 1:numel(spec{m,4})
    bname = spec{m,4}{bb}; e = spec{m,6}(bb); FB = [FC -og.(bname)];
    xP = E.xP{e}; xA = E.xA{e}; PA = [E.oP{e}.P E.oP{e}.A; E.oA{e}.P E.oA{e}.A];
    s = m16_optimum_shift(Xg, FC, FB, p.lb, p.ub, xP, xA, PA);
    % MCDM on each front: equal weights and entropy weights
    DC = -FC(s.maskC,:).*1; DB = -FB(s.maskB,:);          % larger-is-better after sign flip
    tc = m14_mcdm(DC, true(1,size(DC,2)), ones(1,size(DC,2))/size(DC,2));
    tb = m14_mcdm(DB, true(1,size(DB,2)), ones(1,size(DB,2))/size(DB,2));
    te = m14_mcdm(DB, true(1,size(DB,2)), []);
    XC = Xg(s.maskC,:); XB = Xg(s.maskB,:); xc = XC(tc.best,:); xb = XB(tb.best,:); xe = XB(te.best,:);
    Dcomp = norm((xb - xc)./(p.ub - p.lb))/sqrt(numel(p.lb));
    L.shift(end+1,:) = [m, bb, s.nC, s.nB, s.O_P, s.phi_B, s.HV_C, s.HV_B, s.dHV, s.OSI, s.dP, s.dA, Dcomp, tc.stability, tb.stability, te.stability];
    L.shiftlab{end+1} = sprintf('%s / %s', p.name, bname);
    write_csv(fullfile(cfg.out, sprintf('st03_fronts_%s_%s.csv', p.name, bname)), [p.vars(:)' kept {'inC','inB'}], ...
        [Xg cell2mat(cellfun(@(f) og.(f), kept, 'UniformOutput', false)) s.maskC s.maskB]);
    % engineering validation table rows: conventional optimum, benchmark optimum, compromise (conv), compromise (BA)
    D = [xP; xA; xc; xb]; od = p.eval(D); Bv = od.(['B' bname(2:end)]);
    lab = {'conventional optimum (max P)', ['benchmark optimum (max ' bname ')'], 'TOPSIS compromise, conventional front', 'TOPSIS compromise, benchmark-aware front'};
    for q = 1:4
      L.table(end+1,:) = [m bb q D(q,:) nan(1,3-numel(p.lb)) od.P(q) Bv(q) od.P(q)/Bv(q) 1-od.P(q)/Bv(q) od.U(q) od.L(q) od.cv(q)];
      L.tablab{end+1} = sprintf('%s / %s / %s', p.name, bname, lab{q});
    end
  end
end
L.shifthead = {'model','bench','n_Pareto_C','n_Pareto_BA','O_P','phi_B','HV_C','HV_BA','dHV','OSI','dP_rel','dA_rel','D_compromise','stab_C','stab_BA','stab_BA_entropy'};
write_csv(fullfile(cfg.out, 'st03_shift.csv'), L.shifthead, L.shift);
L.tabhead = {'model','bench','row','x1','x2','x3','P','B','A','G','U','L','cv'};
write_csv(fullfile(cfg.out, 'st03_validation_table.csv'), L.tabhead, L.table);
fid = fopen(fullfile(cfg.out, 'st03_labels.txt'), 'w'); fprintf(fid, '%s\n', L.shiftlab{:}); fprintf(fid, '--\n'); fprintf(fid, '%s\n', L.tablab{:}); fclose(fid);
end
