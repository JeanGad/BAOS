function varargout = pc_gpr(action, varargin)
% PC_GPR  Physics-constrained GPR: B is modelled directly, attainment through
% z = logit(A) and losses through log(L), so that 0 < P_hat < B_hat and
% L_hat > 0 hold by construction (constraint satisfaction, not a PINN).
switch action
  case 'fit'
    [X, P, B, L, lb, ub, seed] = varargin{:}; A = min(max(P./B, 1e-6), 1 - 1e-6);
    s.mB = m08_surrogate_models('fit', 'GPR', X, B, lb, ub, seed);
    s.mZ = m08_surrogate_models('fit', 'GPR', X, log(A./(1 - A)), lb, ub, seed);
    s.mL = m08_surrogate_models('fit', 'GPR', X, log(L), lb, ub, seed); varargout{1} = s;
  case 'predict'
    [s, X] = varargin{:}; B = m08_surrogate_models('predict', s.mB, X);
    z = m08_surrogate_models('predict', s.mZ, X); A = 1./(1 + exp(-z));
    varargout{1} = A.*B; varargout{2} = B; varargout{3} = exp(m08_surrogate_models('predict', s.mL, X)); varargout{4} = A;
end
end
