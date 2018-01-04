function gamma_nu  = autocovariance(x, lag_nu)
% auto-covariance
% parameters : lag lag_nu (ndarray) of the samples x = x_1, ..., x_N
% returns :  auto-covariance gamma_nu (ndarray of size lag_nu)

N = length(x);

mhat = mean(x);
cr = (x(1:(N - lag_nu))-mhat) .* (x(1+lag_nu:N) - mhat);
gamma_nu = mean(cr) ;

end