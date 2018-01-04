function X = sample_EulerMaruyama_linearDrift(N, dt, X0)
% Euler-Maruyama discretization
% parameters : int N, number of steps
%              double dt, time step size
%              initial condition X0 of size (d,1)
% return : trajectory X, array (1,N)
    
  % dimesion
  d = length(X0);
  
  % initialize array of samples
  X = zeros(d, N);
    
  % initial condition
  X(1,:) = X0;
    
  % intialize randon numbers, see help randn, help rand
  Z = randn(d, N);
    
    
  for n = 1 : N - 1
                
        X(:,n+1) = X(:,n) + X(:,n) * dt + sqrt(dt) * Z(:,n);
        
  end
    
end