%% Euler-Maruyana discretisation of the SDE
% 
% compute numerical solution of the continuous SDE: dX_t = X_t dt + B_t,
% t>0, X_0=1, where B_t is Wiener process
% Euler-Maruyama discretization: X_n+1 = X_n + dt*X_n + sqrt(dt) * Z_n,
% where Z_n ~ N(0,1) 

% see matlab help for the used functions

%%
% clean the working space
clear all;
close all;
clc;

addpath([pwd,'/functions']);

% choose initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);

myFontSize = 14;

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf(' Euler-Maruyana discretisation of the SDE\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

% intial condition
X0 = 1;

%% Part 1
%% Trajectories

% step size
dt = 0.0001;

% number of trajectories
numberOfTrajectories = 10;

% number of steps
N = 10000;

%lag  = 0:19;

f9 = figure(9);
hold on

for nrtraj = 1 : numberOfTrajectories
    
    
    fprintf('Sampling trajectory = %d\n', nrtraj);
       
    X  = sample_EulerMaruyama_linearDrift(N, dt, X0);
    
    plot((1:N) * dt, X, '-b', 'LineWidth', 2)
    xlabel('time', 'FontSize', myFontSize)
    ylabel('Solution', 'FontSize', myFontSize)
    set(gca, 'FontSize', myFontSize)

end
%%
print(f9,'figures/figure9','-dpng')

%% Tasks:
% 1) Turn off the noise and compare with the exact solution for the
% corresponding ODE

%% Part 2
%% Histogram of X_t at t = 1
% reproduce figure 10 from section 6.3
