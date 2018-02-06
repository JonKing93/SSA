function[] = SSAconvergence(s)
%% Makes a figure to demonstrate MC-SSA convergence (or lack thereof)
%
% -- In --
% s: The output from SSA_Analysis
%

% Get the number of iterations
nIters = size(s.iterSigVals,1);

% Get the iteration data only for the significant singular values
sigIters = s.iterSigVals(:, s.isSigVal);

% Standardize to allow use of same axes
sigIters = zscore(sigIters);

% Plot
figure(); hold on;
plot(1:nIters, sigIters);
title('Convergence of the upper tail of significant singular values.');
xlabel('Number of Monte Carlo Iterations');
ylabel('Standardized Significance');

% Plot the true confidence level of each iteration
figure(); hold on;
plot(1:nIters, s.iterTrueConf);
title('True confidence level for iterative upper tail significance threshold');
xlabel('Number of Monte Carlo Iterations');
ylabel('True confidence level');
