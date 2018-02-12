function[ax] = ssasignificance( s )
%% Makes a semilogy of data singular values vs Monte Carlo values and
% highlights significant values
%
% ----- Inputs -----
%
% s: The output of the SSA_Analysis function
%
% ----- Outputs -----
% 
% ax: An axes handle for the semilogys
%
% ----- Written By -----
%
% Jonathan King, 2018, University of Arizona, jonking93@email.arizona.edu

nVals = numel(s.singVals);

figure;
subplot(2,1,1);
hold on;

% Get legend handles
xval = s.singPeriod(1);
b = semilogy([xval,xval], [s.sigThreshold(1), s.lowerTail(1)], 'k');
dblue = semilogy( xval, s.singVals(1), 'bd');
dred = semilogy( xval, s.singVals(1), 'rd');

% For each singular value
for k = 1:nVals
    
    % Get the x value
    xval = s.singPeriod(k);
    
    % semilogy the 95% bar
    semilogy( [xval, xval], [s.sigThreshold(k), s.lowerTail(k)], 'k');
    
    % Get the color (significance) of each singular value.
    if s.isSigVal(k)
        semilogyStr = 'rd';
    else
        semilogyStr = 'bd';
    end
    
    % semilogy the data singular values
    semilogy( xval, s.singVals(k), semilogyStr);
end

% Labels, legend etc.
xlabel('Singular Vector Period');
ylabel('Singular Value');
title('Comparison of data and Monte Carlo singular values');
legend([b,dblue,dred], sprintf('Monte Carlo %0.f Confidence Interval', 100*(1-s.sig_p)), 'Not Significant', 'Significant');



% REpeat fr frequency
subplot(2,1,2);
hold on;

% Get legend handles
xval = s.singFreq(1);
b = semilogy([xval,xval], [s.sigThreshold(1), s.lowerTail(1)], 'k');
dblue = semilogy( xval, s.singVals(1), 'bd');
dred = semilogy( xval, s.singVals(1), 'rd');

% For each singular value
for k = 1:nVals
    
    % Get the x value
    xval = s.singFreq(k);
    
    % semilogy the 95% bar
    semilogy( [xval, xval], [s.sigThreshold(k), s.lowerTail(k)], 'k');
    
    % Get the color (significance) of each singular value.
    if s.isSigVal(k)
        semilogyStr = 'rd';
    else
        semilogyStr = 'bd';
    end
    
    % semilogy the data singular values
    semilogy( xval, s.singVals(k), semilogyStr);
end

% Labels, legend etc.
xlabel( sprintf('Singular Vector Frequency (cycles / %0.f points)', s.metadata{1,2}) );
ylabel('Log10 Singular Value');
title('Comparison of data and Monte Carlo singular values');
legend([b,dblue,dred], sprintf('Monte Carlo %0.f Confidence Interval', 100*(1-s.sig_p)), 'Not Significant', 'Significant');

end