run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
th = struct('K1',190,'G1',180,'r1',3.21,'K2',3.89,'G2',1.30,'r2',1.20);
lb = [0.05 0.02]; ub = [0.60 0.50];
k = 201; [a,b] = ndgrid(linspace(lb(1),ub(1),k), linspace(lb(2),ub(2),k)); Xg = [a(:) b(:)];
rng(7); X = lb + rand(3000,2).*(ub-lb);
names = {'BT','BA','BS'}; labs = {'B_T = K1/rho1 (constant)','B_A = K_HS+(f)/rho(f)','bound-span attainment'};
fid = fopen('hs_results.csv','w'); fprintf(fid,'bench,class,Rrev_d,Rrev_top,rhoS,Gamma,Gamma_top,P_xP,A_xP,P_xA,A_xA,dA,dP,OSI,f_P,t_P,f_A,t_A,nPareto,thm4_n,thm4_ok\n');
for j = 1:3
  ev = @(Z) addA(model_hs_composite(Z, th), names{j});
  o = ev(X); r = m04_benchmark_independence(o.P, o.B, struct('delta',0.01,'gamma0',0.05), X);
  [xP, oP] = lex_opt(ev, 'P', 'A', lb, ub, Xg, 1e-9); [xA, oA] = lex_opt(ev, 'A', 'P', lb, ub, Xg, 1e-9);
  dA = (oA.A - oP.A)/oA.A; dP = (oP.P - oA.P)/oP.P; osi = norm((xA - xP)./(ub - lb))/sqrt(2);
  og = ev(Xg); m = nd_filter([-og.P -og.A]);
  g = boit_pareto_geometry_check(@(x) ev(x), Xg(m,:), lb, ub);
  printf('%-28s %s Rd=%.3f top=%.3f rho=%.3f G=%.3f/%.3f  xP=(%.3f,%.3f) P=%.2f A=%.3f | xA=(%.3f,%.3f) P=%.2f A=%.3f  dA=%.3f dP=%.3f OSI=%.3f nPar=%d thm4 %d/%d\n', ...
     labs{j}, r.class, r.Rrev_delta, r.Rrev_top, r.rhoS, r.Gamma, r.Gamma_top, xP, oP.P, oP.A, xA, oA.P, oA.A, dA, dP, osi, nnz(m), round(g.frac_ok*g.n), g.n);
  fprintf(fid,'%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%.3f\n', j, r.class, r.Rrev_delta, r.Rrev_top, r.rhoS, r.Gamma, r.Gamma_top, oP.P, oP.A, oA.P, oA.A, dA, dP, osi, xP, xA, nnz(m), g.n, g.frac_ok);
  if j == 2, oPx = model_hs_composite(xP, th); oAx = model_hs_composite(xA, th);
    printf('  gap at xP: context BT-BA=%.2f, microstructure BA-P=%.2f | at xA: %.2f, %.2f\n', oPx.BT-oPx.BA, oPx.BA-oPx.P, oAx.BT-oAx.BA, oAx.BA-oAx.P); end
end
fclose(fid);
% grid data for figure
og = model_hs_composite(Xg, th); dlmwrite('hs_grid.csv', [Xg og.P og.BA og.P./og.BA (og.K-og.KL)./(og.KU-og.KL)], 'precision', 8);
