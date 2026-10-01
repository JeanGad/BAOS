function sel = m06_objective_selection(Y, names, sense, opt)
% M06_OBJECTIVE_SELECTION  Remove objectives that are (practically) strictly
% increasing functions of an already-retained objective (Theorem 1 generalised).
% Y: n x m matrix of candidate objective values over a feasible sample;
% sense: +1 maximise, -1 minimise. Order of columns = physical priority order.
% An objective j is dropped if, oriented for maximisation, it has no delta-resolved
% rank reversal with some retained objective i (R_rev(delta) <= r0).
if nargin < 4, opt = struct(); end
delta = 0.01; r0 = 0; if isfield(opt,'delta'), delta = opt.delta; end
if isfield(opt,'r0'), r0 = opt.r0; end
m = size(Y,2); Z = Y.*sense;                  % orient: larger is better
Zn = (Z - min(Z))./(max(Z) - min(Z) + eps);   % range-normalised (resolution delta on range)
Rd = zeros(m); rho = eye(m); G = zeros(m);
for i = 1:m
  for j = 1:m
    if i == j, continue; end
    Rd(i,j) = resolved_reversal(Zn(:,i), Zn(:,j), delta);
    rho(i,j) = spearman_rho(Z(:,i), Z(:,j));
    [~, o] = sort(Z(:,i)); dz = diff(Z(o,j)); G(i,j) = sum(dz.^2)/(2*(size(Z,1)-1)*var(Z(:,j)) + eps);
  end
end
keep = false(1,m); keep(1) = true; reason = cell(1,m); reason{1} = 'primary objective';
for j = 2:m
  red = find(keep & (Rd(:,j)' <= r0), 1);
  if isempty(red), keep(j) = true; reason{j} = 'independent (reversals with all retained)';
  else, reason{j} = sprintf('redundant: increasing function of %s (Theorem 1)', names{red}); end
end
sel = struct('names',{names},'keep',keep,'reason',{reason},'Rrev_delta',Rd,'rhoS',rho,'Gamma',G);
end
function f = resolved_reversal(a, b, delta)
n = numel(a); if n > 1500, rng(1); k = randperm(n,1500); a = a(k); b = b(k); n = 1500; end
da = a - a.'; db = b - b.'; U = triu(true(n),1);
f = nnz(((da > delta & db < -delta) | (da < -delta & db > delta)) & U)/nnz(U);
end
