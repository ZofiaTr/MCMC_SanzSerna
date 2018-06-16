%% Random Walk Metropolis Algorithm
% The following code implements one dimensional Metropolis random walk, i.e. for $N$ simulation steps, the new state $X^{n+1}$ is obtained from the previous state $X^{n}$ using the proposal $\tilde{X}^{n+1} = X^n + h Z^n,$
% where $h>0$ is the step size, $Z_n$ are independent identically distributed continuous random variables in $\bf{R}^d$ with $d=1$ according to a normal (Gaussian) distribution $\mathcal{N}(0,1)$. The new state is then obtained by Metropolis rule, which defines the acceptance probability $a$ for $\tilde{X}^{n+1}$ as
% $$a  = {\rm min}\left(1, \frac{\rho(\tilde{X}^{n+1})}{\rho(X^n)}\right).$$
% 
% If the proposal is accepted, we set $X^{n+1}=\tilde{X}^{n+1}$, otherwise if the proposal is refused, we set $X^{n+1}=X^{n}$. 
% 
% In the following example, we use this algorithm to sample from the Boltzmann distribution with density
% $$\rho(x) = Z^{-1}{\rm e}^{-\beta x^4},$$
% where $\beta$ is the inverse temperature and $Z$ is the normalization constant such that
% $\int_{\bf{R}} \rho(x)dx = 1$.
% 

% Clean the working space.
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

fprintf('Sampling Random Walk\n')

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
fprintf('Sampling done.\n')

%% 
% Note that the Metropolis rule tests the condition that the proposal is rejected. Since we have already assigned the proposal to the array $X$ at index $n+1$ in the line before, if the proposal is accepted, there is nothing to be done. Testing the refusal condition is more efficient, because in the case of refusal we only assign the value from the previous step to $X_{n+1}$.
%%
% We can plot the generated trajectory $X$ over $N$ steps. Note that in order to make all figures consistent, it is a good idea to define parameters, which, for example, fix the font size. 
myFontSize = 14;

figure(666)
plot(1:length(X),X)
xlabel('steps')
ylabel('X')
%%
% We histogram the samples to see the probability of each state. We save the figure using the print function.

f6 = figure(6);
h = histogram(X, 20, 'Normalization','probability');
xlabel('X', 'FontSize', myFontSize)
ylabel('Probability', 'FontSize', myFontSize)
set(gca, 'FontSize', myFontSize)

print(f6,'figures/figure6','-dpng')
%%
% We fit the target density which is proportional to ${\rm e}^{-x^4}$ in the histogram.
% In order to compute the normalization constant. Since this is only a one
% dimensional problem, we can use quadratures to compute the integral $Z =
% \int_{{\bf{R}}} {\rm e }^{-x^4}dx$.


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

