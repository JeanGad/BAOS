function r = tiedrank1(x)
% TIEDRANK1  Average ranks (toolbox-free replacement for tiedrank).
x = x(:); n = numel(x); [xs, o] = sort(x); r = zeros(n,1); i = 1;
while i <= n
  j = i; while j < n && xs(j+1) == xs(i), j = j + 1; end
  r(o(i:j)) = (i + j)/2; i = j + 1;
end
end
