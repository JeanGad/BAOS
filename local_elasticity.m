function E = local_elasticity(fun, x, h)
% LOCAL_ELASTICITY  d ln Y / d ln x_i by central differences (fun: X -> Y row).
d = numel(x); y0 = fun(x); E = zeros(numel(y0), d);
for i = 1:d
  xp = x; xm = x; xp(i) = x(i)*(1+h); xm(i) = x(i)*(1-h);
  E(:,i) = ((log(fun(xp)) - log(fun(xm)))/(log(1+h) - log(1-h)))';
end
end
