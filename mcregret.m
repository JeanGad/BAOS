addpath(genpath('/home/claude/work_bapimo_0930'));
p = m01_define_problem('brayton'); th0 = p.th; Xg = m07_doe(41^3, p.lb, p.ub, 'grid', 1);
N = 150; rng(11); S = [0.85+0.01*randn(N,1), 0.88+0.01*randn(N,1), 298.15+5*randn(N,1), 0.02+0.02*rand(N,1)];
out = zeros(N,4); bn = {'BT','BA'};
for s = 1:N
  th = th0; th.etac = S(s,1); th.etat = S(s,2); th.T0 = S(s,3); th.dph = S(s,4);
  for j = 1:2
    ev = @(X) mkA(model_brayton(X, th), bn{j});
    [~,oP] = lex_opt(ev,'P','A',p.lb,p.ub,Xg,1e-9); [~,oA] = lex_opt(ev,'A','P',p.lb,p.ub,Xg,1e-9);
    out(s,2*j-1:2*j) = [(oA.A-oP.A)/oA.A, (oP.P-oA.P)/oP.P];
  end
end
dlmwrite('mcregret.csv', out, 'precision', 6);
for k = 1:4, v = sort(out(:,k)); printf('%d median %.4f  CI [%.4f, %.4f]  P(>1%%) %.2f\n', k, median(v), v(round(0.025*N)), v(round(0.975*N)), mean(v > 0.01)); end
