function [X, rejections] = sample_HMC(N, dt, T, V, dV , X0, beta)
% HMC
% parameters : int N, number of steps
%              double dt, time step size
%              potential V, function
%              dV gradient of the potential V, function
%              initial condition X0 of size (d,1)
% return : trajectory X, array (d,N)

rejections =0 ;

% determine dimension from the initial condition
d = length(X0);
% initialize array of samples
X = zeros(d, N);
X(:,1) = X0;

% initialize array of momenta
p = zeros(d, N);

% Hamiltonian defined for X, p \in R^d
H = @(q,p) V(q) + 0.5*p'*p;

U = rand(1, N);

% number of inner steps
L = floor(T/dt);
if (L ==0)
    L=1;
end


for n = 1 : N - 1
    
    %     % sample random momentum
    pRand = normrnd(0,1, size(X0));%randn(d, 1);
    
    xn = X(:,n);
    pn = pRand / beta;
    
    % save previous state
    xnOld = xn;
    pnOld = pn;
    
    dt_rand = dt * rand(1);
    L = floor(T/dt_rand);
    if (L ==0)
        L=1;
    end
    
    for ni = 1: L
        % proposal
        [X_proposal, p_proposal] = sample_Verlet(2, dt_rand, dV, xn, pn);
        xn = X_proposal(:,end);
        pn = p_proposal(:,end);
    end
    
    acceptanceProba =  min(1, exp(-beta*H(xn,pn) + beta*H(xnOld,pnOld)));
    
    if ( acceptanceProba < U(n+1))
        % refuse
        X(:,n+1) = xnOld;
        % reverse momentum
        p(:,n+1) = -pnOld;
        rejections = rejections+1;
    else
        % accept
        X(:,n+1) = xn;
        p(:,n+1) = pn;
        
    end
    
end

kinTemp = mean(p.*p) ;
kinTemp = mean(kinTemp');
fprintf('HMC: kinetic temperature is %f\n', kinTemp );

end