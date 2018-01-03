function X =  sample_Langevin(N, dt, dV, X0, p0)
% Langevin discretization
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
   
  gamma=1;
  kBT=1;
  mass=1;
  a2 = exp(-gamma*dt);
  b2=sqrt(kBT*(1-a2*a2)*mass);

    
  for n = 1 : N - 1
        
        pn = p(:,n);
        xn = X(:,n);
        
        pn = pn - 0.5 * dV(xn) * dt;
        xn = xn + 0.5 * pn * dt;
        
        pn = a2.*pn + b2.*normrnd(0,1, size(pn)) ;
        
        xn = xn + 0.5 * pn * dt;
        pn = pn - 0.5 * dV(xn) * dt;
        
        p(:,n+1) = pn;
        X(:,n+1) = xn;
        
  end
  
  kinTemp = mean(p.*p) ;
kinTemp = mean(kinTemp');
fprintf('Langevin: kinetic temperature is %f\n', kinTemp );
  
end