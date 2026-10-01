function [B, A] = m03_benchmark_operator(out, bname)
% M03_BENCHMARK_OPERATOR  Returns benchmark B and attainment A = P/B for the
% named benchmark field of a model output; checks validity P <= B.
B = out.(bname); A = out.P./B;
if any(out.P > B*(1 + 1e-9)), warning('BAPIMO:validity','P exceeds benchmark %s for %d designs', bname, nnz(out.P > B*(1+1e-9))); end
end
