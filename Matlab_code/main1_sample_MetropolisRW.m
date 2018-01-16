%% Random Walk Metropolis Algorithm
% Proposal by random walk: Y_n+1 = Y_n + Z_n, where Z_n are independent 
% identically distributed continuous random variables in R^d; here d = 1 
% Z_n is a normal (Gaussian) distribution N(m,C) 
% d = 1, hence the density of the normal distirbution is exp(-(x-m)^2)/(2*sigma^2)/(sqrt(2*pi)*sigma), 
% with variance sigma^2 

% see matlab help for the used functions

%% consider N(0,h^2 Id), Z_n ~ N(0,Id)

% clean the working space
clear all;
close all;
clc;


fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('Random Walk Metropolis Algorithm\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');


% number of steps 
N = 1000000;

% step size
h = 1.0;
%inverse  temperature 
beta = 1.0;

% potential
V = @(x) x.^4;

% target probability density
rho = @(x) exp(-beta .*V(x));

% choose initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);

% initialize array of samples
X = zeros(1, N); 

% intialize randon numbers, see help randn, help rand

Z = randn(1, N);
U = rand(1, N);

fprintf('Sampling RW\n')

for n = 1 : N - 1
    
    % proposal 
    X(n+1) = X(n) + h * Z(n);
    
    %Metropolis step: acceptance ratio
    acceptanceRatio = rho(X(n+1)) / rho(X(n));
    
    % random number 
    
    if (acceptanceRatio < U(n+1))
        % refuse
        X(n+1) = X(n);
    end
    
end

fprintf('sampling done\n')

%% show histogram

myFontSize = 14;

f6 = figure(6);
h = histogram(X, 20, 'Normalization','probability');
xlabel('X', 'FontSize', myFontSize)
ylabel('Probability', 'FontSize', myFontSize)
set(gca, 'FontSize', myFontSize)

print(f6,'figures/figure6','-dpng')

%% tasks: 
% 1) fit the target density on the histogram, create plot of X over steps
% 2) comment out the Metropolis step, and rerun the sampling only with the
% random walk proposal, create plot of X over steps

%% Answers:
%% 1) fit the target density

figure(66)
h = histogram(X, 20, 'Normalization','pdf');
hold on
binEdges = h.BinEdges;
% center the bins
binEdges = (binEdges(1:end-1)+binEdges(2:end))/2.0;
normalizationConstant = integral(rho, -inf, inf);
plot(binEdges, rho(binEdges)/normalizationConstant, 'LineWidth', 3)
xlabel('X', 'FontSize', myFontSize)
ylabel('Density', 'FontSize', myFontSize)
set(gca, 'FontSize', myFontSize)

%% plot X over steps
figure(666)
plot(1:length(X),X)
xlabel('steps')
ylabel('X')


