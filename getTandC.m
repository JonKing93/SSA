function[T, C, tsm0] = getTandC(ts, M, algorithm)
%% Builds a Trajectory matrix using either the BK (Broomhead King) or VG (Vautard Ghil) algorithm,
% and constructs the associated covariance matrix.
%
% [T, C] = getTandC(tsm0, M, algorithm)
%
%
% ----- Inputs -----
%
% tsm0: A time series with the means removed. 
%
% M: The size of window used to build SSA Trajectories
%
% algorithm:
%   'BK': Broomhead King
%   'VG': Vautard Ghil 
%
%
% ----- Outputs -----
% 
% T: The Trajectory matrix.
%
% C: The covariance matrix associated with each method

% Run an error check and get pre-allocation sizes
[N] = setup(ts, M);

% Remove mean from time series
tsm0 = detrend(ts, 'constant');

% For the VG algorithm, the window slides off the ends of the time series.
if strcmpi(algorithm, 'VG')   
    
    % Preallocate
    ncols = N+M-1;
    T = NaN(M, ncols);
    C = NaN(M,M);

    % Build the trajectory matrix
    for k = 1:M-1    % Before the window completely covers the series
        T(M-k+1:M, k) = tsm0(1:k);
    end
    for k = M: ncols-M   % When there is complete overlap
        T(:,k) = tsm0(k-M+1:k);
    end
    for k = ncols-M+1:ncols  % When the window has fallen off the series
        T(1:ncols-k+1,k) = tsm0(k-M+1:N);
    end

    % Calculate covariance matrix using Toeplitz symmetry.
    for k = 1:M
        notNAN = ( ~isnan(T(1,:)) & ~isnan(T(k,:)) ); 
        C(1,k) = sum( T(1,notNAN) .* T(k,notNAN)) ./ sum(notNAN); % The first row / column
    end
    C(:,:) = toeplitz( C(1,:));   

% For the BK algorithm, the window stops at the edges of the time series.
elseif strcmpi(algorithm, 'BK')
 
    % Preallocate
    ncols = N-M+1;
    T = NaN(M, ncols);
    C = NaN(M,M);

    % Fill matrix
    for k = 1:ncols  % For each value in the window M
        T(:,k) = tsm0(k:M+k-1);
    end

    % Calculate covariance matrix
    C(:,:) = cov(T');
        
else
    error('Unrecognized algorithm');
end
    
    
end

% ----- Helper Functions -----
function[N] = setup(tsm0, M)

% Check that tsm0 is a vector
if ~isvector(tsm0) || any(isnan(tsm0(:)))
    error('Data_m0 must be a vector and cannot contain NaN');
end

% Get size of Data
N = length(tsm0);

% Check that the window size is positive
if M < 0
    error('M must be positive');
elseif mod(M,1) ~= 0
    error('M must be an integer');
end

% Check that the window size does not exceed the time series size
if M > N
    error('M exceeds time series length');
end

end




 