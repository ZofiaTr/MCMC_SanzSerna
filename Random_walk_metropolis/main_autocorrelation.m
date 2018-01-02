%% Autocorrelation of the Random Walk Metropolis Algorithm
% 
% compute the empirical auto-covariance 

% see matlab help for the used functions

%%
% clean the working space
clear all;
close all;
clc;


% choose initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);


myFontSize = 14;

% step size
h = [0.5, 1, 2, 4];

% number of various stepsizes
nrh = length(h);

lag  = 1:20;

figure(7)
hold on

i=0;

for stepSize = h
    
    i = i+1;

    [X, rejections] = sample_MetropolisRW( stepSize);
    
    rho_nu = compute_empirical_auto_correlation_coeff(X, lag);
    
    subplot(2,2,i)
    plot(lag, rho_nu, 'LineWidth', 2)
    xlabel('lag', 'FontSize', myFontSize)
    ylabel('Correlation', 'FontSize', myFontSize)
    ylim([0,1])
    title(['h = ', num2str(stepSize), ', ', num2str(rejections), ' rejections'])
    set(gca, 'FontSize', myFontSize)

end

%% Tasks:
% Compare with matlab's autocorr function.

%% Solution
%  [rho_nu, lag, ~] = autocorr(X, 20);








