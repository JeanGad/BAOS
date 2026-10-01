run('/home/claude/BAOS_toolbox_rel/setup_baos.m'); addpath('/home/claude/rev');
for n = [100 300 1000 3000 10000]
  rng(1); X = rand(n,2); [P,B] = known_case(X,4);
  tic; for k=1:3, r = m04_benchmark_independence(P,B,struct('delta',0.01),X); end; t=toc/3;
  printf('n=%d  t=%.3f s  class %s\n', n, t, r.class);
end
