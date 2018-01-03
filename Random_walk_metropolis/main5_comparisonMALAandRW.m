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
fprintf('Compare MALA and Random-walk\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');

% choose initial seed, comment out to turn off, see help rng
% seed=0;
% rng(seed);
myFontSize = 14;

lag  = 0:19;

% define the target density :
k = 100.0;
d = 100;
% radius
r = @(x)  norm(x, 2);
% potential: notation from the paper- mathcal{L} = -V
V = @(x) 0.5 * k * (r(x) - 1)^2; %debug 0.5* x'*x;%
% gradient of the potential
%debug with dV = @(x) x;%
dV = @(x) k * (r(x) - 1).* x/r(x);%k * (r(x) - 1).*sign(x);
% target probability density
rho = @(x) exp(- V(x));

% initial condition 
X0 = ones(d, 1);

%% Random Walk
% step size
h = 0.025;
numberOfSteps = 10000;

f11 = figure(11);
hold on

fprintf('Sampling with RW\n');
 
[X, rejections] = sample_MetropolisRW(numberOfSteps, h, rho , X0);
    
rho_nu = compute_empirical_auto_correlation_coeff(X(1,:), lag);

subplot(1,2,1)
scatter(X(1,:), X(2,:),  2)
xlabel('X_1', 'FontSize', myFontSize)
ylabel('X_2', 'FontSize', myFontSize)
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])


subplot(1,2,2)
plot(lag, rho_nu, 'LineWidth', 2)
xlabel('lag', 'FontSize', myFontSize)
ylabel('Correlation', 'FontSize', myFontSize)
ylim([0,1])
xlim([0,max(lag)])
title(['h = ', num2str(h),', ',num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
set(gca, 'FontSize', myFontSize)

%%
virial = mean(X.*dV(X)) ;
confTempRW = mean(virial');
fprintf('RW: configurational temperature is %f\n', confTempRW );

%%
print(f11,'figures/figure11','-dpng')

%%
% %% MALA
% % step size

h = 0.0000001;
numberOfSteps = 10000;

fprintf('Sampling with MALA\n');
 
[X, rejections] = sample_MALA(numberOfSteps, h, V, dV , X0);
    
rho_nu = compute_empirical_auto_correlation_coeff(X(1,:), lag);

%%
virial = mean(X.*dV(X)) ;
confTempMALA = mean(virial');
fprintf('MALA: configurational temperature is %f\n', confTempMALA );

%%
f12 = figure(12);
hold on

subplot(1,2,1)
scatter(X(1,:), X(2,:),  2)
xlabel('X_1', 'FontSize', myFontSize)
ylabel('X_2', 'FontSize', myFontSize)
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])


subplot(1,2,2)
plot(lag, rho_nu, 'LineWidth', 2)
xlabel('lag', 'FontSize', myFontSize)
ylabel('Correlation', 'FontSize', myFontSize)
%ylim([0,1])
xlim([0,max(lag)])
title(['h = ', num2str(h),', ',num2str(rejections), ' rejections, rate ', num2str(rejections / length(X))])
set(gca, 'FontSize', myFontSize)


%%
print(f12,'figures/figure12','-dpng')

