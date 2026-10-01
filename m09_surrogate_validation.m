function m = m09_surrogate_validation(y, yh)
% M09_SURROGATE_VALIDATION  R2, RMSE, MAE, MAPE (%) of predictions yh vs truth y.
e = yh - y; m.R2 = 1 - sum(e.^2)/sum((y - mean(y)).^2); m.RMSE = sqrt(mean(e.^2));
m.MAE = mean(abs(e)); m.MAPE = 100*mean(abs(e)./max(abs(y), eps)); m.maxAPE = 100*max(abs(e)./max(abs(y), eps));
end
