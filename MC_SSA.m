function[surrVals] = MC_SSA(ts, M, algorithm, singVecs, MC, noise, varargin)
%% Runs a Monte Carlo singular spectrum analysis. Generates surrogate eigenvalues.
%
% [surrVals, iterSigVals, iterTrueConf] = MC_SSA(tsm0, M, algorithm, singVecs, MC, noise)
% Conducts an MC-SSA significance test on a set of singular vectors.
% Returns the set of surrogate singular values as well as surrogate
% eigenvalues at the significance level for each iteration.
%
% MC_SSA(..., 'parallel')
% To speed up runtime, runs the Monte Carlo process in parallel using 
% MATLAB's default settings. If the Distributed Computing Toolbox is not 
% licensed, reverts to serial computing. 
%
% MC_SSA(..., 'showProgress')
% Displays a progress bar showing the progress through the Monte Carlo 
% iterations. This option is disabled when computing in parallel.
%
% MC_SSA(..., 'estimateRuntime')
% Estimates total runtime based on the number of iterations and available
% workers in the computing pool.
%
%
% ----- Inputs -----
% 
% ts: A time series.
%
% M: The embedding dimension / sampling window.
%
% algorithm: The algorithm used to build SSA trajectories
%       'BK': Broomhead-King
%       'VG': Vautard-Ghil
%
% singVecs: The singular vectors for the time series.
%
% MC: The number of Monte Carlo iterations in the MC-SSA
% 
% noise: A flag for the noise used in the Rule N significance test
%       'red' (Default): AR(1) noise with added white noise
%       'white': white Gaussian noise. (mean = 0, variance = 1)
%
%
% ----- Outputs -----
%
% surrVals: The matrix of surrogate singular values

% Parse inputs. Error Check. Setup.
[parallel, showProgress, estimateRuntime] = setup(ts, MC, varargin{:});

% Preallocate surrogate singular values
surrVals = NaN(MC, M);

% Initialize the parallel pool if computing in parallel
nWorkers = 1;
if parallel
    pool = gcp;
    nWorkers = pool.NumWorkers;
end

% Initialize the progress bar if displaying
if showProgress
    progressbar(0);
end

% Run the first iteration. Estimate runtime if desired
if estimateRuntime
    startTime = tic;
end
surrVals(1,:) = mcssaStep( ts, M, algorithm, singVecs, noise );
if estimateRuntime
    time = toc(startTime);
    time = time*MC/nWorkers;
    h = time / 360;
    m = mod(time,360)/60;
    s = mod(time,60);
    fprintf('Estimated runtime: %0.f hour(s), %0.f minute(s), %0.f seconds \r\n',h,m,s);
end

% Run the MC-SSA
if parallel             % In parallel
    parfor k = 2:MC
        surrVals(k,:) = mcssaStep( ts, M, algorithm, singVecs, noise );
    end
    
else                    % In serial
    for k = 2:MC
        surrVals(k,:) = mcssaStep( ts, M, algorithm, singVecs, noise );
        
        % Update progress bar if displaying
        if showProgress
            progressbar(k/MC);
        end
    end
end

end


%%%%% Helper functions %%%%%%
function[parallel, showProgress, estimateRuntime] = setup(ts, MC, varargin)

% Parse the optional inputs
[parallel, showProgress, estimateRuntime] = parseInputs( varargin, ...
{'parallel','showProgress','estimateRuntime'},{false, false, false},{'b','b','b'} );

% Ensure tsm0 is a vector
if ~isvector(ts)
    error('ts must be a vector');
end

% Ensure MC is positive
if MC < 1 || mod(MC,1)~=0
    error('Monte Carlo number must be a positive integer');
end

% If parallel, check for Distributed Computing Toolbox
if parallel
    if ~license('test', 'Distrib_Computing_Toolbox')
        warning( sprintf('The Distributed Computing Toolbox is not licensed on this computer.\r\nCannot run in parallel. Reverting to default serial mode...')); %#ok<SPWRN>
        parallel = false;
    end
end

% Disable 'showProgress' for parallel computing
if parallel && showProgress
    warning('Cannot display progress for parallel computations.');
    showProgress = false;
end

end

function[surrVals] = mcssaStep(ts, M, algorithm, singVecs, noise)

% Build a surrogate time series
surr = randNoiseSeries(noise, ts, 1);

% Get the covariance matrix
[~,surrC] = getTandC( surr, M, algorithm);

% Project onto data eigenvector basis
surrProj = singVecs' * surrC * singVecs;

% Extract the diagonals elements. These are the surrogate singular values
surrVals = diag( surrProj );

end

