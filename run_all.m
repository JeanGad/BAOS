% RUN_ALL  BAPIMO: complete, reproducible pipeline (MATLAB R2025b / GNU Octave 8+).
% No toolboxes required. Results -> results/*.csv|.mat ; figures -> figures/ via
% m19_plot_results (MATLAB) or plot_results.py (used for the reported figures).
% Total run time on one CPU core in Octave 8.4: roughly 15-20 minutes.
clear; close all;
root = fileparts(mfilename('fullpath')); addpath(genpath(root));
cfg = struct('seed', 20260930, ...     % master seed; every stage derives its seeds from it
             'nBOIT', 3000, ...         % feasible random sample size for BOIT
             'delta', 0.01, ...         % relative resolution (1%) used in BOIT and selection
             'gamma0', 0.05, ...        % functional-dependence threshold for Gamma
             'nTrain', 120, ...         % surrogate training size (maximin LHS)
             'surrogate', 'GPR', ...    % chosen by the rule in st04 (see doc, Part V)
             'nRuns', 30, 'nSobol', 4096, 'nMC', 2000, ...
             'out', fullfile(root, 'results'));
if ~exist(cfg.out, 'dir'), mkdir(cfg.out); end
S  = st01_synthetic(cfg);                 % Level 1-2 validation: known-answer suite
E  = st02_engineering_boit(cfg);          % BOIT + exact optima + Theorem 4 / Cor. 4.1 check
save('-v7', fullfile(cfg.out, 'st02.mat'), 'E');
L  = st03_levels_shift(cfg, E);           % objective selection, Models 1-3, optimum shift, MCDM
R  = st04_surrogates(cfg);                % surrogate comparison (accuracy, Q2, extrapolation, physics)
V  = st04b_surrogate_moo(cfg, E);         % m15: surrogate MOO + high-fidelity re-evaluation + UQ
A2 = st05_algorithms(cfg, 2);             % NSGA-II / MOEA/D / MOPSO, 30 runs, statistics
A3 = st05_algorithms(cfg, 3);
Sb = st06_sensitivity(cfg, E);            % Sobol + Proposition 5 check
Rb = st07_robustness(cfg, E);             % Monte Carlo UQ, delta method, chance-constrained robust MOO
Ab = st08_ablation(cfg);                  % objective-set and benchmark-definition ablation
st09_export(cfg, E);                      % figure data
if exist('OCTAVE_VERSION', 'builtin') == 0, m19_plot_results(cfg); end
disp('BAPIMO pipeline complete.');
