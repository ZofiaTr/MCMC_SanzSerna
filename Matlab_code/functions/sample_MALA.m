function [X, rejections] = sample_MALA(N, h, V, dV , X0)
% MALA
% parameters : int N, number of steps
%              double h, time step size
%              potential V, function
%              dV gradient of the potential V, function
%              initial condition X0 of size (d,1)
% return : trajectory X, array (d,N)

dt = 0.5*h^2;
fprintf('Mala step size is %f\n', dt)
rejections = 0 ;

% determine dimension from the initial condition
d = length(X0);
% initialize array of samples
X = zeros(d, N); 
X(:,1) = X0;

%U = rand(1, N);

for n = 1 : N - 1
    
    xn = X(:,n);
    xnOld = xn;
    
    % proposal 
    Z = normrnd(0,1, size(xn));
    
    xn = xn - 0.5 * dV(xn) * dt + sqrt(dt) * Z;
    
    X(:,n+1) = xn;
    
    T_new =  V(xn) ;
    T_old =  V(xnOld) ;
      
    acceptanceProbaX =  exp( - (T_old - T_new) ) * ( TproposalMALA(xn,xnOld, dV ,dt) / TproposalMALA(xnOld, xn, dV, dt));
    
    %Metropolis step: acceptance ratio
    
    if (  acceptanceProbaX < rand(1) )
        % refuse
        X(:,n+1) = xnOld;
        rejections = rejections+1;
 
       
    end


    
end


     

end