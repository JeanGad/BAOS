run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
o = struct('delta',0.01,'gamma0',0.05,'nboot',200,'seed',3);
fid = fopen('twotier.csv','w'); fprintf(fid,'case,global,local,local_ci,Gamma,G_lo,G_hi,Gamma_top,Gt_lo,Gt_hi,Rrev_top,Rrev,R_lo,R_hi\n');
pr = @(lab,r) fprintf(fid,'%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',lab,r.class,r.class_top,r.class_top_ci,r.Gamma,r.Gamma_ci,r.Gamma_top,r.Gamma_top_ci,r.Rrev_top,r.Rrev_delta,r.Rrev_ci);
rng(20260930); X = rand(3000,2);
for c = 1:6, [P,B] = known_case(X,c); r = m04_benchmark_independence(P,B,o,X); pr(sprintf('S%d',c), r); end
labs = {'Brayton B_T','Brayton B_C','Brayton T3 fixed B_T','Brayton T3 fixed B_C','Heat exchanger','Solar tau-alpha','Solar F_R tau-alpha'};
for c = 1:7
  M = csvread(sprintf('/home/claude/work_bapimo_0930/results/st02_sample_%d.csv', c), 1, 0);
  P = M(:,end-2); B = M(:,end-1); r = m04_benchmark_independence(P,B,o,M(:,1:end-3)); pr(labs{c}, r);
end
th = struct('K1',190,'G1',180,'r1',3.21,'K2',3.89,'G2',1.30,'r2',1.20); lb = [0.05 0.02]; ub = [0.60 0.50];
rng(7); X = lb + rand(3000,2).*(ub-lb); m = model_hs_composite(X, th);
nm = {'BT','BA','BS'}; lc = {'Composite B_T','Composite B_HS','Composite span'};
for j = 1:3, r = m04_benchmark_independence(m.P, m.(nm{j}), o, X); pr(lc{j}, r); end
fclose(fid); disp('done')
