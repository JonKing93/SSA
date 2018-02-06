function[RCs, signals] = getRCs(singVecs, T, algorithm)
%% Calculates the RCs for a SSA output
% 
% [RCs, signals] = getRCs(singVecs, T, algorithm)
%
%
% ----- Inputs -----
%
% singVecs: The singular vectors for a trajectory matrix.
%
% T: A trajectory matrix.
%
% algorithm: The algorithm used to calculate the trajectory matrices
%   'BK': Broomhead - King
%   'VG': Vautard - Ghil
%
%
% ----- Outputs -----
% 
% RCs: The reconstructed components for a trajectory matrix and associated
%       singular values. Each column is an RC.
%
% signals: The signals calculated and used to construct the RCs. Each
%       column is a signal

% Error check, get some initial sizes
[M, N] = setup(singVecs, T, algorithm);

% Preallocate
RCs = NaN(M, N);

% Get the signals
signals = singVecs' * T;
    
% Get the RCs by convolving the singular vectors with the signals
for k=1:M
    notNAN = ~isnan(signals(k,:));
    RCs(k,:) = conv( singVecs(:,k), signals(k, notNAN) );
end

% Apply convolution weights to normalize the RCs
for k=1:N
    % Case 1: incomplete overlap of functions because we are at the front
    % of the vector
    if k <= M-1
        RCs(:,k) = RCs(:,k) ./ k;

    % Case 2: Complete overlap of the two functions
    elseif k >= M && k <= N-M+1
        RCs(:,k) = RCs(:,k) ./ M;

    % Case 3: Incomplete overlap because we are at the end of the vector
    else
        RCs(:,k) = RCs(:,k) ./ (N-k+1);
    end
end

% Return with individual RCs in columns
RCs = permute(RCs, [2,1]);
signals = permute(signals, [2,1]);

end

% ----- Helper functions -----
function[M, N] = setup(singVecs, T, algorithm)

% Ensure inputs are matrices
if ~ismatrix(singVecs) || ~ismatrix(T) || hasNaN(singVecs) || hasNaN(T)
    error('singVecs and T must be matrices and cannot contain NaN');
end

% Calculate some sizes
[M, ncols] = size(T);

% Get the number of points
switch algorithm
    case 'VG'
        N = ncols - M + 1;
    case 'BK'
        N = ncols + M - 1;
    otherwise
        error('Unrecognized algorithm');
end

end