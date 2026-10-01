function A = st05_algorithms(cfg, which, algsel)
% ST05_ALGORITHMS  NSGA-II vs MOEA/D vs MOPSO on the high-fidelity Brayton
% problem, 30 independent runs, equal budget (100 x 100 evaluations).
% which = 2: max[P, A_A];  which = 3: max[P, A_A, U].
p = m01_define_problem('brayton');
if which == 2, nm = {'P','AA'}; else, nm = {'P','AA','U'}; end
fun = @(X) m02_physical_model(X, p, nm, ones(1,numel(nm))); M = numel(nm);
algs = {@nsga2, @moead, @mopso}; an = {'NSGA-II','MOEA/D','MOPSO'}; R = cfg.nRuns;
fr = cell(R,3); tm = zeros(R,3); trc = cell(R,3);
% provisional reference from dense grid for the HV trace (M = 2 only)
Xg = m07_doe(41^3, p.lb, p.ub, 'grid', 1); [Fg, cg] = fun(Xg); Fg = Fg(cg <= 0,:); Fg = Fg(nd_filter(Fg),:);
id0 = min(Fg,[],1); nd0 = max(Fg,[],1);
if nargin < 3, algsel = 1:3; end
cache = fullfile(cfg.out, sprintf('st05_cache_M%d.mat', M));
if exist(cache, 'file'), load(cache); end
for a = algsel
  for r = 1:R
    if ~isempty(fr{r,a}), continue; end
    opt = struct('pop', 100, 'gen', 100, 'seed', cfg.seed + 1000*a + r, 'M', M);
    if M == 2, opt.trace = @(F) hv((F - id0)./(nd0 - id0), [1.1 1.1]); end
    tic; res = feval(algs{a}, fun, p.lb, p.ub, opt); tm(r,a) = toc;
    fr{r,a} = res.F; if M == 2, trc{r,a} = res.trace; end
    save('-v7', cache, 'fr', 'tm', 'trc');
  end
  printf('  st05(M=%d): %s done\n', M, an{a}); fflush(stdout);
  save('-v7', cache, 'fr', 'tm', 'trc');
end
if any(cellfun(@isempty, fr(:))), A = []; return; end
Fall = [Fg; cell2mat(fr(:))]; Fref = Fall(nd_filter(Fall),:); id = min(Fref,[],1); nd = max(Fref,[],1);
HV = zeros(R,3); IGD = HV; GD = HV; SP = HV; NS = HV;
for a = 1:3, for r = 1:R
  q = m11_pareto_analysis(fr{r,a}, Fref, id, nd); HV(r,a) = q.HV; IGD(r,a) = q.IGD; GD(r,a) = q.GD; SP(r,a) = q.SP; NS(r,a) = size(fr{r,a},1);
end, end
A.an = an; A.HV = HV; A.IGD = IGD; A.GD = GD; A.SP = SP; A.T = tm; A.NS = NS; A.Fref = Fref;
A.stHV = m18_statistics(HV, an); A.stIGD = m18_statistics(IGD, an); A.stT = m18_statistics(tm, an);
tag = sprintf('M%d', M);
write_csv(fullfile(cfg.out, ['st05_runs_' tag '.csv']), {'HV_NSGA2','HV_MOEAD','HV_MOPSO','IGD_NSGA2','IGD_MOEAD','IGD_MOPSO', ...
   'GD_NSGA2','GD_MOEAD','GD_MOPSO','SP_NSGA2','SP_MOEAD','SP_MOPSO','t_NSGA2','t_MOEAD','t_MOPSO'}, [HV IGD GD SP tm]);
st = [A.stHV A.stIGD A.stT]; S = zeros(numel(st), 4);
for i = 1:numel(st), S(i,:) = [ceil(i/3) st(i).p st(i).padj st(i).A12]; end
write_csv(fullfile(cfg.out, ['st05_tests_' tag '.csv']), {'metric(1=HV,2=IGD,3=time)','p','p_holm','A12'}, S);
fid = fopen(fullfile(cfg.out, ['st05_tests_' tag '_pairs.txt']), 'w');
for i = 1:numel(st), fprintf(fid, '%s vs %s\n', st(i).a, st(i).b); end; fclose(fid);
write_csv(fullfile(cfg.out, ['st05_reffront_' tag '.csv']), nm, -Fref);
if M == 2
  T = zeros(100, 3); for a = 1:3, T(:,a) = mean(cell2mat(trc(:,a)'), 2); end
  write_csv(fullfile(cfg.out, 'st05_hvtrace_M2.csv'), an, T);
  write_csv(fullfile(cfg.out, 'st05_example_fronts_M2.csv'), {'alg','P','AA'}, ...
     [ones(size(fr{1,1},1),1) -fr{1,1}; 2*ones(size(fr{1,2},1),1) -fr{1,2}; 3*ones(size(fr{1,3},1),1) -fr{1,3}]);
end
save('-v7', fullfile(cfg.out, ['st05_' tag '.mat']), 'A');
end
