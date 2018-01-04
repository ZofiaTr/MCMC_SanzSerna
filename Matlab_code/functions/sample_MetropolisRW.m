function [X, rejections] = sample_MetropolisRW(N, h, rho, X0)
% parameters :
%  
%   number of steps N
%   stepsize h
%   rho is target density
%   initial condition X_0, size(1, d)
% return : trajectory and number of rejections

% determine dimension from the initial condition
d = length(X0);
% initialize array of samples
X = zeros(d, N); 
X(:,1) = X0;
rejections = 0;

% intialize randon numbers, see help randn, help rand

Z = randn(d, N);
U = rand(1, N);

for n = 1 : N - 1
    
    % proposal 
    X(:,n+1) = X(:,n) + h * Z(:,n);
    
    %Metropolis step: acceptance ratio
    acceptanceRatio = rho(X(:,n+1)) / rho(X(:,n));
    
    % random number 
    
    if (acceptanceRatio < U(n+1))
        % refuse
        X(:,n+1) = X(:,n);
        rejections = rejections+1;
    end
    
end


end