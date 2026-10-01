function m19_plot_results(cfg)
% M19_PLOT_RESULTS  MATLAB rendering of the main figures from results/*.csv.
% NOTE: written for MATLAB R2025b; NOT executed in this study (headless Octave
% plotting failed in the build environment). The reported figures were rendered
% by plot_results.py from the same CSV files.
out = cfg.out; fdir = fullfile(fileparts(out), 'figures'); if ~exist(fdir,'dir'), mkdir(fdir); end
S = readmatrix(fullfile(out, 'st09_surface_brayton.csv')); nr = 111; nt = 121;
g = @(c) reshape(S(:,c), nr, nt); rp = g(1); T3 = g(2);
names = {'P','B_T','B_A','A_T','A_A'}; cols = [3 4 5 6 7];
f = figure('Visible','off','Position',[100 100 1200 650]);
for k = 1:5
  subplot(2,3,k); contourf(rp, T3, g(cols(k)), 20, 'LineColor','none'); colorbar; hold on;
  contour(rp, T3, g(10), [1000 1000], 'w', 'LineWidth', 1); title(names{k}); xlabel('r_p'); ylabel('T_3 (K)');
end
subplot(2,3,6); contourf(rp, T3, 1 - g(7), 20, 'LineColor','none'); colorbar; title('G_A');
exportgraphics(f, fullfile(fdir, 'M01_brayton_surfaces.png'), 'Resolution', 200); close(f);
for c = [1 2 5]
  F = readmatrix(fullfile(out, sprintf('st02_front_%d.csv', c))); D = readmatrix(fullfile(out, sprintf('st02_sample_%d.csv', c)));
  f = figure('Visible','off'); scatter(D(:,end-2), D(:,end), 3, [.7 .7 .7]); hold on; scatter(F(:,end-1), F(:,end), 12, 'r', 'filled');
  xlabel('P'); ylabel('A'); title(sprintf('Case %d: Pareto set of max[P,A]', c));
  exportgraphics(f, fullfile(fdir, sprintf('M02_front_case%d.png', c)), 'Resolution', 200); close(f);
end
end
