function [X, rejections] = sample_MALA(N, h, V, dV , X0)
% MALA
% parameters : int N, number of steps
%              double h, time step size
%              potential V, function
%              dV gradient of the potential V, function
%              initial condition X0 of size (d,1)
% return : trajectory X, array (d,N)

dt = sqrt(h);
rejections =0 ;

% determine dimension from the initial condition
d = length(X0);
% initialize array of samples
X = zeros(d, N); 
X(:,1) = X0;

% target probability density
rho = @(x) exp(- V(x));

% dV
dV2 = @(x)0.5*dV(x);

U = rand(1, N);

for n = 1 : N - 1
    
    % proposal 
    X_proposal = sample_EulerMaruyama(2, dt, dV2, X(:,n));
    X(:,n+1) = X_proposal(:,end);
    
    T_new = V( X(:,n+1));
    T_old = V( X(:,n));
    
    alpha = U(n+1);
    
    acceptanceProbaX =  exp( - (T_old - T_new) ) * ( TproposalMALA(X(:,n+1),X(:,n), dV ,dt) / TproposalMALA(X(:,n),X(:,n+1), dV, dt));
    
    %Metropolis step: acceptance ratio
    %acceptanceProbaX = rho(X(:,n+1)) / rho(X(:,n));
            
    if (alpha > acceptanceProbaX)
        % refuse
        X(:,n+1) = X(:,n);
        rejections = rejections+1;
    end
    
end


     

end