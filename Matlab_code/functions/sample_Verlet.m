function [X, p] = sample_Verlet(N, dt, dV, X0, p0)
% Verlet discretization
% parameters : int N, number of steps
%              double dt, time step size
%              dV gradient of the potential V, function
%              initial condition X0 of size (d,1)
%              initial condition for momenta p0 of size (d,1)
% return : 
%              trajectory X, array (d,N)
%              trajectory of momenta p, array (d,N)

  % dimesion
  d = length(X0);
  
  % initialize array of samples
  X = zeros(d, N);
  p = zeros(d, N);
    
  % initial condition
  X(:,1) = X0;
  p(:,1) = p0;
   
    
  for n = 1 : N - 1
        
        pn = p(:,n);
        xn = X(:,n);
        
        pn = pn - 0.5 * dV(xn) * dt;
        xn = xn + pn * dt;
        pn = pn - 0.5 * dV(xn) * dt;
        
        p(:,n+1) = pn;
        X(:,n+1) = xn;
        
  end
    
end