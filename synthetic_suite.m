function [P, B, info] = synthetic_suite(X, c)
% SYNTHETIC_SUITE  Six controlled problems on x in [0,1]^2 with analytically
% known P-B relationships (expected BOIT class and optimum shift in info).
x1 = X(:,1); x2 = X(:,2);
switch c
  case 1 % S1: constant benchmark -> exact redundancy
    P = 1 + x1 + 0.5*sin(pi*x2); B = 3*ones(size(x1));
    info = struct('name','S1 constant B','expect','R0','shift',0);
  case 2 % S2: B = 2 sqrt(P) (elasticity 1/2 < 1) -> monotone redundancy
    P = 0.2 + 2*x1.*(0.5 + 0.5*x2); B = 2*sqrt(P);
    info = struct('name','S2 B=f(P), elasticity<1','expect','R1','shift',0);
  case 3 % S3: B = P/h(P), h hump-shaped -> non-monotone functional dependence
    P = 0.2 + x1 + 0.3*x2; Pn = (P - 0.2)/1.3; h = 0.2 + 0.7*4*Pn.*(1 - Pn);
    B = P./h;
    info = struct('name','S3 B=f(P), non-monotone h','expect','D2','shift',1);
  case 4 % S4: design-dependent gap -> independence with optimum shift
    P = 1 + x1 + 0.6*x2; B = P + 0.2 + 1.2*x2.^2;
    info = struct('name','S4 design-dependent B','expect','I3','shift',1);
  case 5 % S5: benchmark varies by 0.5% only -> redundant at resolution delta = 1%
    P = 1 + x1 + 0.5*x2; B = 3*(1 + 0.005*sin(2*pi*x2));
    info = struct('name','S5 weakly varying B','expect','R1','shift',0);
  case 6 % S6: A independent of P but maxima coincide -> independent, no shift
    P = 0.1 + 0.9*x1; A = 0.3 + 0.6*x2; B = P./A;
    info = struct('name','S6 independent, non-conflicting','expect','I3','shift',0);
end
end
