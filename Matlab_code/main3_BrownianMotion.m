%% Standard Brownian motion
% Random walk: B_n+1 = B_n + sqrt(dt) * Z_n, where Z_n are independent
% identically distributed continuous random variables in R^d; here d = 1
% Z_n is a normal (Gaussian) distribution N(m,C)

% see matlab help for the used functions

%% consider N(0,h^2 Id), Z_n ~ N(0,Id)

% clean the working space
clear all;
close all;
clc;


fprintf('@@@@@@@@@@@@@@@\n');
fprintf('Brownian motion\n');
fprintf('@@@@@@@@@@@@@@@\n');


% number of steps
N = 10000;

% step size
dt = 0.0001;

% choose initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);

% number of trajectories to be simulated
numberOfTrajectories = 2;

myFontSize = 14;
f8 = figure(8);


for nrtraj = 1: numberOfTrajectories
    
    % initialize array of samples
    X = zeros(1, N);
    
    % intialize randon numbers, see help randn, help rand
    Z = randn(1, N);
    
    
    for n = 1 : N - 1
        
        % Brownian motion
        X(n+1) = X(n) + sqrt(dt) * Z(n);
        
    end
    
    plot((1:N) .* dt, X)
    hold on
    xlabel('time', 'FontSize', myFontSize)
    ylabel('B(t)', 'FontSize', myFontSize)
    set(gca, 'FontSize', myFontSize)
    
    
end

%%
print(f8,'figures/figure8','-dpng')

%%