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
       
    X  = sample_EulerMaruyama(N, dt);
    
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

%% Solution
% 1) exact solution of the ODE is exp(t)

%% Part 2
%% Histogram of X_t at t = 1

% number of trajectories
numberOfTrajectories = 100000;

% step size
dt = 0.001;

% number of steps, fix final time t = 1
N = floor(1 / dt);

XfinalTime = zeros(1,numberOfTrajectories);



for nrtraj = 1 : numberOfTrajectories
    
    
    fprintf('Sampling trajectory number %d\n', nrtraj);
       
    X  = sample_EulerMaruyama(N, dt);
    XfinalTime(nrtraj) = X(end);
        
end
%%

f10 = figure(10);
h = histogram(XfinalTime, 15, 'Normalization','pdf');
hold on

mu = exp(1);
sigma = sqrt((exp(1)^2 - 1)/2);

plot(XfinalTime,normpdf(XfinalTime,mu,sigma),'*', 'MarkerSize', 1)
xlabel('X', 'FontSize', myFontSize)
ylabel('Frequencies', 'FontSize', myFontSize)
set(gca, 'FontSize', myFontSize)

print(f10,'figures/figure10','-dpng')
