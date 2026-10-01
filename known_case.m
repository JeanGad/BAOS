function [P, B] = known_case(X, c)
x1 = X(:,1); x2 = X(:,2);
switch c
  case 1, P = 1 + x1 + 0.5*sin(pi*x2); B = 3*ones(size(x1));
  case 2, P = 0.2 + 2*x1.*(0.5 + 0.5*x2); B = 2*sqrt(P);
  case 3, P = 0.2 + x1 + 0.3*x2; Pn = (P - 0.2)/1.3; B = P./(0.2 + 2.8*Pn.*(1 - Pn));
  case 4, P = 1 + x1 + 0.6*x2; B = P + 0.2 + 1.2*x2.^2;
  case 5, P = 1 + x1 + 0.5*x2; B = 3*(1 + 0.005*sin(2*pi*x2));
  case 6, P = 0.1 + 0.9*x1; B = P./(0.3 + 0.6*x2);
end
end
