function X = sample_EulerMaruyama(N, dt, dV, X0)
% Euler-Maruyama discretization
% parameters : int N, number of steps
%              double dt, time step size
%              dV gradient of the potential V, function
%              initial condition X0 of size (d,1)
% return : trajectory X, array (d,N)
    
  % dimesion
  d = length(X0);
  
  % initialize array of samples
  X = zeros(d, N);
    
  % initial condition
  X(:,1) = X0;
    
  % intialize randon numbers, see help randn, help rand
  Z = randn(d, N);
    
    
  for n = 1 : N - 1
        % attention: recall the notation mathcal{L} = -V
        X(:,n+1) = X(:,n) - dV(X(:,n)) * dt + sqrt(dt) * Z(:,n);
        
  end
    
end