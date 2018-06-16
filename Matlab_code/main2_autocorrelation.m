%% Autocorrelation of the Random Walk Metropolis Algorithm
% We compute the empirical auto-covariance 
% Implement auto-covariance function which takes array of samples $X$ and
% lag array and returns empirical auto-covariance computed as
% $$=\mu_{\nu} = \frac{\gamma_{\nu}}{\gamma_0}$ where 
% $\gamma_{\nu} = \frac{1}{N+1}\sum_{i=0}^{N-\nu}(x_i - \hat{m})(x_{i+\nu} - \hat{m}), \quad \hat{m} = \frac{1}{N+1}\sum_{i=0}^N x_i$.


%% 
% Before starting the main code, we wrap-up the sampling part of the Metropolis Random-Walk, that we have seen in the previous section, into a separate function and save it into the folder "functions". This function takes number of steps to be performed $N$, steps size $h$, distribution that should be sampled $\mu$ and initial state $X_0$ and returns an array of $N$ samples $X$ and an array of rejections of size $N$.

% clean the working space
clear all;
close all;
clc;

addpath([pwd,'/functions']);

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('Autocorrelation of the Random Walk Metropolis Algorithm\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

% fix initial seed, comment out to turn off, see help rng
seed=0;
rng(seed);
myFontSize = 14;

%inverse  temperature 
beta = 1.0;
% potential
V = @(x) x.^4;
% target probability density
mu = @(x) exp(-beta .*V(x));

%% 
% We want to investigate several step sizes, we therefore have an array of step sizes $h$, which will be iterated over later on.
% step size
h = [0.5, 1, 2, 4];
% number of various stepsizes
nrh = length(h);
% number of steps 
N = 1000000;
%initial condition 
X0=0;
%%
% In order to compute the empirical auto-covariance, we choose the lag array. We also initialize a figure to allow writing inside during the iterations over the step sizes.
lag  = 0:19;

f7 = figure(7);
hold on

%%
% What follows is the main for-loop over the step sizes: for every step size $h$, we perform the sampling and compute the empirical auto-correlation. The implementation of this function compute_empirical_auto_correlation_coeff will be shown at the end of this section.
% Finally, within each iteration, we plot the auto-covariance in a subplot belonging to the pre-intialized figure.
% We also compare our auto-correlation function with matlab's function
% autocorr and plot the results in one figure.
% intialize a counter
i=0;
for stepSize = h
    
    fprintf('Sampling with step size h = %f\n', stepSize);
    i = i+1;

    [X, rejections] = sample_MetropolisRW(N, stepSize, mu, X0);
    %%%%% compute_empirical_auto_correlation_coeff
    rho_nu = compute_empirical_auto_correlation_coeff(X, lag);
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % MATLAB's autocorrelation function
    [rho_nu_matlab, lag_matlab, ~] = autocorr(X, 20);
    
    subplot(2,2,i)
    plot(lag, rho_nu, '-b', 'LineWidth', 2)
    hold on
    plot(lag_matlab, rho_nu_matlab, '--r', 'LineWidth', 2)
    xlabel('lag', 'FontSize', myFontSize)
    ylabel('Correlation', 'FontSize', myFontSize)
    ylim([0,1])
    xlim([0,max(lag)])
    legend({'Our','MATLAB'})
    title(['h = ', num2str(stepSize), ', rejection rate ', num2str(rejections / length(X))])
    set(gca, 'FontSize', myFontSize)

end
% save figure
print(f7,'figures/figure7','-dpng')

%%
% We have implemented the empirical auto-covariance using the two following functions:

% 
% <include>functions/compute_empirical_auto_correlation_coeff.m</include>
%

% 
% <include>functions/autocovariance.m</include>
%