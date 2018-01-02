function gamma_nu  = autocovariance(x, lag_nu)

N = length(x);

mhat = mean(x);
gamma_nu = sum( (x(1:(N - lag_nu))-mhat) .* x(1+lag_nu:N) - mhat) / (N+1);

end