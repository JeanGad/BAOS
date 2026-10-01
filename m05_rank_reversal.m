function r = m05_rank_reversal(P, A, delta, nmax, seed)
% M05_RANK_REVERSAL  Pairwise rank-reversal statistics between P and A.
% R_rev        : fraction of pairs with P_i>P_j and A_i<A_j (discordant pairs)
% R_rev_delta  : delta-resolved reversals (log-separation > delta in BOTH P and A)
% tau          : Kendall tau-a.  Pairs are over a random subsample of size <= nmax.
if nargin < 4, nmax = 1500; end
if nargin < 5, seed = 1; end
P = P(:); A = A(:); n = numel(P);
if n > nmax, rng(seed); idx = randperm(n, nmax); P = P(idx); A = A(idx); n = nmax; end
dP = P - P.'; dA = A - A.';
U = triu(true(n),1); Npairs = nnz(U);
conc = (dP.*dA > 0) & U; disc = (dP.*dA < 0) & U;
lP = log(P); lA = log(A); dlP = lP - lP.'; dlA = lA - lA.';
discd = ((dlP > delta & dlA < -delta) | (dlP < -delta & dlA > delta)) & U;
r.Rrev = nnz(disc)/Npairs; r.Rrev_delta = nnz(discd)/Npairs;
r.tau = (nnz(conc) - nnz(disc))/Npairs; r.Npairs = Npairs; r.n = n;
end
