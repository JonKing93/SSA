function[ax] = ssaconvergence(s)
%% Makes a figure to demonstrate MC-SSA convergence (or lack thereof)
%
% [ax] = ssaconvergence(s)
%
% ----- Inputs -----
%
% s: The output from SSA_Analysis
%
% ----- Outputs -----
%
% ax: A handle to axes for each of the plots
%
% ----- Written By -----
%
% Jonathan King, 2018, University of Arizona, jonking93@email.arizona.edu

if ~isfield(s,'MC_p') || ~isfield(s,'MCsigThresh') || ~isfield(s,'MCisSig')
    warning('Insufficient data for ssaconvergence. Try running a convergence test...');
    return;
end

% Significance level
figure();
plot( s.MC_p);
xlabel('Number of Monte Carlo Iterations');
ylabel('True signficance level of signifcance test.');
title('True Significance Level Tested at each Monte Carlo Iteration');
ax = gca;

% Singular values
figure()
plot( zscore( s.MCsigThresh(:,s.isSigVal) ) );
xlabel('Number of Monte Carlo Iterations');
ylabel('Normalized Significance Threshold of Significant Singular Values');
title('Significance Threshold of Significant Singular Values for Monte Carlo Iterations');
ax = [ax;gca];

% Significance Pass/Fail
figure()
imagesc( s.MCisSig );
cmap = [[1,1,1];[0,0,1]];
colormap(cmap);
pfalse = patch(0,0,[1,1,1]);
ptrue = patch(0,0,[0,0,1]);
legend([ptrue, pfalse], 'Significant', 'Not Significant');
xlabel('Singular Value');
ylabel('Monte Carlo Iteration');
title('Significant Singular Values at each Monte Carlo Iteration');
ax = [ax; gca];
end




% % Get the number of iterations
% nIters = size(s.iterSigVals,1);
% 
% % Get the iteration data only for the significant singular values
% sigIters = s.iterSigVals(:, s.isSigVal);
% 
% % Standardize to allow use of same axes
% sigIters = zscore(sigIters);
% 
% % Plot
% figure(); hold on;
% plot(1:nIters, sigIters);
% title('Convergence of the upper tail of significant singular values.');
% xlabel('Number of Monte Carlo Iterations');
% ylabel('Standardized Significance');
% 
% % Plot the true confidence level of each iteration
% figure(); hold on;
% plot(1:nIters, s.iterTrueConf);
% title('True confidence level for iterative upper tail significance threshold');
% xlabel('Number of Monte Carlo Iterations');
% ylabel('True confidence level');
% 
% end