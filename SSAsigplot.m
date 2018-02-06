function[s] = SSAsigplot( s, xType, periodUnit)
%% Makes a plot of data eigenvalues vs Monte Carlo values and plots the
% desired significance line. Determines the significant frequencies and
% periods.
%
% [sigFreq, sigPeriod] = SSA_Sig( s )
%
% ----- Inputs -----
%
% s: The output of the SSA_Analysis function
%
% xType: Specifies the variables for the x-axis
%   'Freq': Plots frequency along the x-axis
%   'Period': Plots period along the x-axis
%
% periodUnit: The number of observations in a period unit 
%       (e.g. periodUnit = 12 for monthly observations on an annual scale)
%
% ----- Outputs -----
%
% s: The following fields will be added to the input structure
%   ei   
%   sigFreq: The frequencies of siginificant oscillations
%
%   sigPeriod: The periods of significant oscillations
%
% ---
% Jonathan King, 2017

nseries = size(s.singVals,2);

if strcmpi(xType, 'freq')
    xVar = s.maxFreq;
    xStr = 'Frequency';
elseif strcmpi( xType, 'period')
    xVar = s.maxPeriod / periodUnit;
    xStr = 'Period';
else
    error('Unrecognized xType');
end

figure();clf;hold on;

% Plot the range of surrogate eigenvalues within the confidence
% interval tails of the Monte Carlo tests
for m = 1:size(s.singVals,1)
    ax1 = semilogy([xVar(m), xVar(m)], [s.lowSigVals(m), s.upSigVals(m)], 'k');
end

% Add the data eigenvalues
ax2 = semilogy( xVar, s.singVals, 'rd');

% Make the significant eigenvalues blue
ax3 = semilogy( xVar(s.isSigVal), s.singVals(s.isSigVal), 'bd');

% Add labels etc.
legend([ax1,ax2,ax3],'Monte Carlo 95% Confidence Interval','Data Eigenvalues','Significant Eigenvalues');
title('Comparison of data eigenvalues with monte carlo eigenvalues.');
xlabel(sprintf('Eigenvalue %s', xStr));
ylabel('Eigenvalue Power');

end