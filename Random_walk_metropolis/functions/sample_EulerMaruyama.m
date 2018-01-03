function X = sample_EulerMaruyama(N, dt)
% Euler-Maruyama discretization
% parameters : int N, number of steps
%              double dt, time step size
% return : trajectory X, array (1,N)

  % initialize array of samples
    X = zeros(1, N);
    
    % initial condition
    X(1) = 1.0;
    
    % intialize randon numbers, see help randn, help rand
    Z = randn(1, N);
    
    
    for n = 1 : N - 1
                
        X(n+1) = X(n) + X(n) * dt + sqrt(dt) * Z(n);
        
    end
    
end