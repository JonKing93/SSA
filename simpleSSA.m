function[singVals, singVecs, tsm0, T, C] = simpleSSA(ts, M, algorithm)
%% Runs a single-spectrum analysis.
%
% [singVals, singVecs, ts_m0, T, C] = simpleSSA(Data, M)
% performs singular spectrum analysis on each time series in a set of data
% series.
%
% [...] = simpleSSA(..., algorithm)
% performs singular spectrum analysis using a trajectory matrix constructed
% via a particular algorithm.
%
%
% ----- Inputs -----
%
% ts: A time series vector with equally spaced observations.
%
% M: The embedding dimension.
%
% algorithm:
%       'BK': Broomhead-King -- Slightly less bias for nonstationary time series
%       'VG': Vautard-Ghil -- Enhanced noise reduction for short time series
%
%
% ----- Outputs -----
%
% singVals: The singular values of each SSA of each time series
%
% singVecs: The singular vectors of each SSA of each time series 
%
% ts_m0: a time series with the mean removed
%
% T: The trajectory matrix for the time series.
%
% C: The covariance matrix for the trajectory matrix.

% Get trajectory matrix and associated covariance
[T, C, tsm0] = getTandC(ts, M, algorithm);

% Run an svd of the covariance matrix
[~,S, singVecs] = svd(C);

% Get the singular values from the eigenvalues of the C matrix
singVals = diag(S);

% Make the singular vectors majority positive elements. (This is optional).
for k = 1:size(singVecs,2)
    if sum( singVecs(:,k)<0 ) > (size(singVecs,1)/2)
        singVecs(:,k) = -1 * singVecs(:,k);
    end
end

end