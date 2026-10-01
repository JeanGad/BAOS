function rho = spearman_rho(a, b)
% SPEARMAN_RHO  Spearman rank correlation with average ranks for ties.
ra = tiedrank1(a(:)); rb = tiedrank1(b(:));
ra = ra - mean(ra); rb = rb - mean(rb);
rho = (ra'*rb)/sqrt((ra'*ra)*(rb'*rb) + eps);
end
