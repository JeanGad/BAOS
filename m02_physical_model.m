function [F, cv, out] = m02_physical_model(X, prob, names, sense)
% M02_PHYSICAL_MODEL  Generic objective wrapper: evaluates the replaceable
% physics model and returns the requested responses in minimisation form.
out = prob.eval(X); F = zeros(size(X,1), numel(names));
for j = 1:numel(names), F(:,j) = -sense(j)*out.(names{j}); end
cv = out.cv;
end
