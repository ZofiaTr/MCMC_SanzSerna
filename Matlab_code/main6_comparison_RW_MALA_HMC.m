%% Compare RW, MALA and HMC
% d = 3, k =100
% sample exp(-0.5*k(r-1)^2), r=|x|, x \in R^d
% see matlab help for the used functions

%%
% clean the working space
clear all;
close all;
clc;

addpath([pwd,'/functions']);

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('Compare HMC, MALA and Random-walk\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

% choose initial seed, comment out to turn off, see help rng
% seed=0;
%rng(seed);

myFontSize = 14;

lag  = 0:19;

% define the target density :
k = 100;
d = 3;
% radius
r = @(x) norm(x, 2);
% potential: notation from the paper- mathcal{L} = -V
V = @(x) 0.5 * k * (r(x) - 1)^2; 
% gradient of the potential
diffV = @(x) dV(x, k, d);
% target probability density
rho = @(x) exp(- V(x));

% initial condition 
X0 = zeros(d, 1);
X0(1) = 1; 
%% Random Walk
% step size
h = 0.5;
numberOfSteps = 2000;

f14 = figure(14);
hold on

fprintf('Sampling with RW\n');
 
[X, rejections] = sample_MetropolisRW(numberOfSteps, h, rho , X0);
    
rho_nu = compute_empirical_auto_correlation_coeff(X(1,:), lag);

subplot(1,2,1)
scatter(X(1,:), X(2,:),  2)
xlabel('X_1', 'FontSize', myFontSize)
ylabel('X_2', 'FontSize', myFontSize)
xlim([-2, 2])
ylim([-2, 2])


subplot(1,2,2)
plot(lag, rho_nu, 'LineWidth', 2)
xlabel('lag', 'FontSize', myFontSize)
ylabel('Correlation', 'FontSize', myFontSize)
ylim([0,1])
xlim([0,max(lag)])
title(['RW: h = ', num2str(h),', ',num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
set(gca, 'FontSize', myFontSize)

%%
confTempRW = compute_configurational_temperature(X, diffV);
fprintf('RW: configurational temperature is %f\n', confTempRW );

%%
print(f14,'figures/figure14','-dpng')

%%
% %% MALA
% % step size
h = 0.25;
numberOfSteps = 2000;

fprintf('Sampling with MALA\n');
 
[X, rejections] = sample_MALA(numberOfSteps, h, V, diffV , X0);
    
rho_nu = compute_empirical_auto_correlation_coeff(X(1,:), lag);

%%
confTempMALA = compute_configurational_temperature(X, diffV);
fprintf('MALA: configurational temperature is %f\n', confTempMALA );

%%
f15 = figure(15);
hold on

subplot(1,2,1)
scatter(X(1,:), X(2,:),  2)
xlabel('X_1', 'FontSize', myFontSize)
ylabel('X_2', 'FontSize', myFontSize)
xlim([-2, 2])
ylim([-2, 2])


subplot(1,2,2)
plot(lag, rho_nu, 'LineWidth', 2)
xlabel('lag', 'FontSize', myFontSize)
ylabel('Correlation', 'FontSize', myFontSize)
%ylim([0,1])
xlim([0,max(lag)])
title(['MALA: h = ', num2str(h),', ',num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
set(gca, 'FontSize', myFontSize)


%%
print(f15,'figures/figure15','-dpng')

%% HMC
% step size
dt = 0.1;
T = 1;
numberOfSteps = 1000;%%200;

fprintf('Sampling with HMC\n');

%X  = sample_Langevin(numberOfSteps, dt, diffV , X0, zeros(size(X0)));
[X, rejections] = sample_HMC(numberOfSteps, dt, T, V, diffV , X0);

rho_nu = compute_empirical_auto_correlation_coeff(X(1,:), lag);

%%


confTempHMC = compute_configurational_temperature(X, diffV);
fprintf('HMC: configurational temperature is %f\n', confTempHMC );

%%
f16 = figure(16);
hold on

subplot(1,2,1)
scatter(X(1,:), X(2,:),  2)
xlabel('X_1', 'FontSize', myFontSize)
ylabel('X_2', 'FontSize', myFontSize)
xlim([-2, 2])
ylim([-2, 2])

subplot(1,2,2)
plot(lag, rho_nu, 'LineWidth', 2)
xlabel('lag', 'FontSize', myFontSize)
ylabel('Correlation', 'FontSize', myFontSize)
%ylim([0,1])
%xlim([0,max(lag)])
title(['HMC: T = ', num2str(T),', dt = ', num2str(dt),', ',num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
set(gca, 'FontSize', myFontSize)

%%
print(f16,'figures/figure16','-dpng')
