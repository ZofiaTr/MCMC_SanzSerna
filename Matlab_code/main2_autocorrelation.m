%% Autocorrelation of the Random Walk Metropolis Algorithm
% 
% compute the empirical auto-covariance 
% see matlab help for the used functions

%%
% clean the working space
clear all;
close all;
clc;

addpath([pwd,'/functions']);

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('Autocorrelation of the Random Walk Metropolis Algorithm\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

% choose initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);
myFontSize = 14;

%inverse  temperature 
beta = 1.0;

% potential
V = @(x) x.^4;

% target probability density
rho = @(x) exp(-beta .*V(x));

% step size
h = [0.5, 1, 2, 4];
% number of various stepsizes
nrh = length(h);
% number of steps 
N = 1000000;
%initial condition 
X0=0;

lag  = 0:19;

f7 = figure(7);
hold on

i=0;

for stepSize = h
    
    fprintf('Sampling with step size h = %f\n', stepSize);
    i = i+1;

    [X, rejections] = sample_MetropolisRW(N, stepSize, rho, X0);
    %%%%% compute_empirical_auto_correlation_coeff
    rho_nu = compute_empirical_auto_correlation_coeff(X, lag);
    %%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,2,i)
    plot(lag, rho_nu, 'LineWidth', 2)
    xlabel('lag', 'FontSize', myFontSize)
    ylabel('Correlation', 'FontSize', myFontSize)
    ylim([0,1])
    xlim([0,max(lag)])
    legend(['h = ', num2str(stepSize)])
    title([num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
    set(gca, 'FontSize', myFontSize)

end
%%
print(f7,'figures/figure7','-dpng')

%% Tasks:
% 1) Compare with matlab's autocorr function.
% 2) Try more values of the step size h.

%% Solution
% [rho_nu_matlab, lag_matlab, ~] = autocorr(X, 20);
%     
%     subplot(2,2,i)
%     plot(lag, rho_nu, 'LineWidth', 2)
%     hold on
%     plot(lag_matlab, rho_nu_matlab, '-r','LineWidth', 2)
% 







