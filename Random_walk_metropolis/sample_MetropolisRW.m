function [X, rejections] = sample_MetropolisRW(stepSize)

% number of steps 
N = 1000000;

% step size
h = stepSize;

%inverse  temperature 
beta = 1.0;

% potential
V = @(x) x.^4;

% target probability density
rho = @(x) exp(-beta .*V(x));

% initialize array of samples
X = zeros(1, N); 
rejections = 0;

% intialize randon numbers, see help randn, help rand

Z = randn(1, N);
U = rand(1, N);


for n = 1 : N - 1
    
    % proposal 
    X(n+1) = X(n) + h * Z(n);
    
    %Metropolis step: acceptance ratio
    acceptanceRatio = rho(X(n+1)) / rho(X(n));
    
    % random number 
    
    if (acceptanceRatio < U(n+1))
        % refuse
        X(n+1) = X(n);
        rejections = rejections+1;
    end
    
end


end