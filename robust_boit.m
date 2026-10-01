run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
grp = @(c) find(strcmp(c, {'R0','R1','D2','I3'}));       % 1,2 redundant; 3 D2; 4 I3
exg = [1 1 3 4 1 4]; exs = [0 0 1 1 0 0]; R = 20;
% Study 1: sample size, noise-free
ns = [50 100 200 500 1000 3000]; S1 = zeros(numel(ns), 6, 2);
for a = 1:numel(ns), for c = 1:6, ok = 0; oks = 0;
  for r = 1:R
    rng(1000*a + 10*c + r); X = rand(ns(a),2); [P,B] = known_case(X,c);
    q = m04_benchmark_independence(P,B,struct('delta',0.01,'gamma0',0.05),X);
    g = grp(q.class); g(g==2) = 1; ok = ok + (g == exg(c)); oks = oks + (q.shift == exs(c));
  end
  S1(a,c,:) = [ok oks]/R;
end, end
fid = fopen('study1.csv','w'); fprintf(fid,'n,case,class_acc,shift_acc\n');
for a=1:numel(ns), for c=1:6, fprintf(fid,'%d,%d,%.2f,%.2f\n',ns(a),c,S1(a,c,1),S1(a,c,2)); end, end; fclose(fid);
% Study 2: evaluation noise at n = 1000, fixed delta vs noise-matched delta
sg = [0 0.0025 0.005 0.01 0.02]; S2 = zeros(numel(sg), 6, 4);
for a = 1:numel(sg), for c = 1:6, acc = zeros(1,4);
  for r = 1:R
    rng(5000*a + 10*c + r); X = rand(1000,2); [P,B] = known_case(X,c);
    P = P.*exp(sg(a)*randn(size(P))); B = B.*exp(sg(a)*randn(size(B)));
    for k = 1:2
      if k == 1, dl = 0.01; else, dl = max(0.01, 3*sg(a)); end
      q = m04_benchmark_independence(P,B,struct('delta',dl,'gamma0',0.05),X);
      g = grp(q.class); g(g==2) = 1;
      acc(2*k-1) = acc(2*k-1) + (g == exg(c)); acc(2*k) = acc(2*k) + (q.shift == exs(c));
    end
  end
  S2(a,c,:) = acc/R;
end, end
fid = fopen('study2.csv','w'); fprintf(fid,'sigma,case,class_acc_d1,shift_acc_d1,class_acc_d3s,shift_acc_d3s\n');
for a=1:numel(sg), for c=1:6, fprintf(fid,'%.4f,%d,%.2f,%.2f,%.2f,%.2f\n',sg(a),c,S2(a,c,:)); end, end; fclose(fid);
disp('done')
