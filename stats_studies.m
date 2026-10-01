run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
grp = @(c) find(strcmp(c, {'R0','R1','D2','I3'}));
% ---------- (a) graded-independence family, d = 4 (not used for tuning) ----------
lam = [0 0.02 0.05 0.1 0.2 0.5 1]; R = 20; fa = fopen('graded.csv','w');
fprintf(fa,'lambda,u_pop,sigma,rho,frac_I3,frac_R1,frac_D2\n');
gfam = @(X,l) deal(1 + 1.5*mean(X,2), 0.35 + 0.35*mean(X,2) + l*0.15*(X(:,1)-X(:,2)+X(:,3)-X(:,4))/2);
for l = lam
  rng(1); Xb = rand(2e5,4); [Pb, Ab] = gfam(Xb, l); [~,o] = sort(Pb); u = sum(diff(Ab(o)).^2)/(2*(numel(Ab)-1)*var(Ab));
  for cfg = [0 0; 0.01 0; 0.01 0.8]'
    sg = cfg(1); rh = cfg(2); cnt = zeros(1,4);
    for r = 1:R
      rng(100*r + round(1000*l)); X = rand(1000,4); [P,A] = gfam(X,l); B = P./A;
      e1 = randn(size(P)); e2 = rh*e1 + sqrt(1-rh^2)*randn(size(P));
      P = P.*exp(sg*e1); B = B.*exp(sg*e2);
      q = m04_benchmark_independence(P,B,struct('delta',max(0.01,3*sg),'gamma0',0.05,'rtol',0.005*(sg>0)),X);
      cnt(grp(q.class)) = cnt(grp(q.class)) + 1;
    end
    fprintf(fa,'%.2f,%.4f,%.3f,%.1f,%.2f,%.2f,%.2f\n', l, u, sg, rh, cnt(4)/R, (cnt(1)+cnt(2))/R, cnt(3)/R);
  end
end
fclose(fa);
% ---------- (b) correlated noise on S1-S6 ----------
exg = [1 1 3 4 1 4]; fb = fopen('corrnoise.csv','w'); fprintf(fb,'rho,rule,acc\n');
for rh = [0 0.5 0.9], for rule = 1:2
  acc = 0;
  for c = 1:6, for r = 1:R
    rng(7000 + 10*c + r); X = rand(1000,2); [P,B] = known_case(X,c);
    e1 = randn(size(P)); e2 = rh*e1 + sqrt(1-rh^2)*randn(size(P)); P = P.*exp(0.01*e1); B = B.*exp(0.01*e2);
    if rule == 1, op = struct('delta',0.01,'gamma0',0.05); else, op = struct('delta',0.03,'gamma0',0.05,'rtol',0.005); end
    q = m04_benchmark_independence(P,B,op,X); g = grp(q.class); g(g==2) = 1; acc = acc + (g == exg(c));
  end, end
  fprintf(fb,'%.1f,%d,%.3f\n', rh, rule, acc/(6*R));
end, end
fclose(fb);
% ---------- (c) high-dimensional problems, d = 10 and 20 ----------
fc = fopen('highdim.csv','w'); fprintf(fc,'d,case,class_acc,dA_true,dA_found,dP_true,dP_found,xd_true,xd_found\n');
for d = [10 20]
  P_ = @(X) 1 + mean(X(:,1:d-1),2) + 0.6*X(:,d);
  cases = {@(X) 3*ones(size(X,1),1)+P_(X)*0, @(X) 2*sqrt(P_(X)), @(X) P_(X) + 0.2 + 1.2*X(:,d).^2};
  ex = [1 1 4];
  for c = 1:3
    ok = 0;
    for r = 1:R
      rng(9000 + 100*d + 10*c + r); X = rand(1000,d); P = P_(X); B = cases{c}(X);
      if c == 1, B = 3.2*ones(size(P)); end
      q = m04_benchmark_independence(P,B,struct('delta',0.01,'gamma0',0.05),X); g = grp(q.class); g(g==2) = 1; ok = ok + (g == ex(c));
    end
    if c == 3
      t = linspace(0,1,100001)'; Aa = (2+0.6*t)./(2.2+0.6*t+1.2*t.^2); [Am,i] = max(Aa); xdt = t(i);
      PA = 2 + 0.6*xdt; dAt = (Am - (2.6/(2.6+1.4)))/Am; dPt = (2.6 - PA)/2.6;
      ev = @(Z) struct('P',P_(Z),'A',P_(Z)./cases{3}(Z),'cv',zeros(size(Z,1),1));
      rng(5); Xs = rand(2000,d); lbv = zeros(1,d); ubv = ones(1,d);
      [xP,oP] = lex_opt(ev,'P','A',lbv,ubv,Xs,1e-9); [xA,oA] = lex_opt(ev,'A','P',lbv,ubv,Xs,1e-9);
      dAf = (oA.A - oP.A)/oA.A; dPf = (oP.P - oA.P)/oP.P; xdf = xA(d);
    else, dAt = 0; dAf = 0; dPt = 0; dPf = 0; xdt = NaN; xdf = NaN; end
    fprintf(fc,'%d,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', d, c, ok/R, dAt, dAf, dPt, dPf, xdt, xdf);
  end
end
fclose(fc); disp('done')
