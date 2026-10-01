run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
grp = @(c) find(strcmp(c, {'R0','R1','D2','I3'})); exg = [1 1 3 4 1 4]; R = 20;
sg = [0 0.0025 0.005 0.01 0.02]; rt = [0.005 0.01]; S = zeros(numel(sg),6,2);
for a = 1:numel(sg), for c = 1:6, acc = [0 0];
  for r = 1:R
    rng(5000*a + 10*c + r); X = rand(1000,2); [P,B] = known_case(X,c);
    P = P.*exp(sg(a)*randn(size(P))); B = B.*exp(sg(a)*randn(size(B)));
    for k = 1:2
      q = m04_benchmark_independence(P,B,struct('delta',max(0.01,3*sg(a)),'gamma0',0.05,'rtol',rt(k)),X);
      g = grp(q.class); g(g==2) = 1; acc(k) = acc(k) + (g == exg(c));
    end
  end
  S(a,c,:) = acc/R;
end, end
fid = fopen('study2b.csv','w'); fprintf(fid,'sigma,case,acc_rtol005,acc_rtol01\n');
for a=1:numel(sg), for c=1:6, fprintf(fid,'%.4f,%d,%.2f,%.2f\n',sg(a),c,S(a,c,:)); end, end; fclose(fid);
% engineering-like check: classes of the Brayton cases are unchanged by rtol? (reported R_rev(delta) >= 0.008)
disp('done')
