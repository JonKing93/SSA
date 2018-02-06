function[s] = SSA_Analysis(ts, M, varargin)
%% Performs singular spectra analyses for a collection of time series.
% Computes RCs (reconstructed components), and performs a Monte Carlo
% significance test
%
% [s] = SSA_Analysis(ts, M)
% Conducts a singular spectrum analysis of a time series using the 
% Broomhead-King algorithm. Applies a 1000 iteration Monte Carlo SSA (MC-SSA) 
% procedure to test singular vector significance (p<0.05) against a red
% AR(1)-noise background. Calculates reconstructed components. Estimates
% the period of singular vectors. Returns a structure 's' with analysis 
% output and metadata. Plots analysis data.
%
% [s] = SSA_Analysis(..., 'noplot')
% Suppresses the output plots.
%
% ANALYSIS OPTIONS:
%
% SSA_Analysis(..., 'algorithm', algorithm)
% Specifies whether to use the Broomhead-King or Vautard-Ghil algorithm.
%
% SSA_Analysis(..., 'MC', MC)
% Specifies the number of Monte Carlo iterations to use during significance
% testing. A greater number of iterations increases statistical robustness,
% but also raises runtime.
%
% SSA_Analysis(..., 'noiseType', noise)
% Specifies wheter to perform the significance test against a red AR(1), or
% white noise process. By default, SSA_Analysis uses a red significance
% test. Use a white significance test (Gaussian noise, variance=1, mean=0)
% for processes with minimal autocorrelation. Use a red test for processes
% with high autocorrelation or long memory. A red test is generally more
% stringent and applicable to many geophysical / geological processes.
%
% SSA_Analysis(..., 'p', p)
% Specify the significance level cutoff to use for significance testing. By
% default SSA_Analysis uses p<0.05 to assign significance.
% 
% SSA_Analysis(..., 'noSigTest')
% Blocks the MC-SSA significance test.
%
% RUNTIME OPTIONS:
% Huge time series? MC-SSA taking forever? Try these options...
%
% SSA_Analysis(..., 'parallel')
% If possible, runs the MC-SSA significance tests in parallel. This can
% significantly speed runtime for calculations longer than several minutes.
% If the "Distributed Computing Toolbox" is not licensed, notifies the user
% and proceeds in serial.
%
% SSA_Analysis(..., 'estimateRuntime')
% Displays an estimate of total runtime for the MC-SSA.
%
% SSA_Analysis(..., 'showProgress')
% Displays a progress bar with percent completion and estimated remaining
% time for MC-SSA. Not available for parallel computations.
%
% SSA_Analysis(..., noConvergeTest')
% A flag to block the test for Monte Carlo convergence. This may slightly
% improve runtime for very large significance tests.
%
%
% ----- Inputs -----
%
% ts: A time series vector with equally spaced observations.
%
% M: The embedding dimension / sampling window length. M corresponds to the 
%       longest singular vector wavelength that can be extracted from the 
%       data. Larger values of M thus provide information on more
%       wavelengths, while smaller values of M yield higher confidence in 
%       extracted information.
%
% algorithm: The algorithm used to compute the SSA
%       'BK': Broomhead-King - Slightly less bias for nonstationary time series
%       'VG': Vautard-Ghil - Enhanced noise reduction for short time series
%
% MC: The number of Monte Carlo iterations used in the Rule N significance
%       test. Must be a positive integer.
%
% noise: A flag for the noise used in the Rule N significance test
%       'red' (Default): AR(1) noise with added white noise
%       'white': white Gaussian noise. (mean = 0, variance = 1)
%
% p: The significance level that the significance test should pass. Must be
%    a positive number on the interval (0, 1).
%
%
% ----- Outputs -----
%
% s: a structure containing some or all of the following fields.
%
%   tsm0: The analysis time series. This is the time series normalized to a
%        mean of 0.
%
%   singVals: The singular values of the time series. 
%
%   singVecs: The singular vectors of each time series. Each column is a
%             singular vector.
%
%   singFreq: An estimate of the frequency of each singular vector.
%
%   singPeriod: An estimate of the period of each singular vector.
%
%   T: The trajectory matrix constructed for the analyses. Each dim1 x
%       dim2 matrix is the trajectory matrix for one time series.
%
%   C: The covariance matrix for the trajectory matrix.
%
%   RCs: The reconstructed components of each time series. Each dim1 x dim2
%       matrix corresponds to a particular time series. Each column of each
%       matrix contains 1 reconstructed component.
%
%   surrVals: The surrogate eigenvalues from the Monte Carlo significance test.
%
%   isSigVal: A boolean/logical vector indicating whether singular values
%   passed the significance test.
%
%   sigThreshold: The threshold singular value for significance. Singular
%       values must be greater than or equal to this threshold to be significant.
%
%   sig_p: The actual significance level of the MC-SSA. (Depending on
%   the number of Monte Carlo iterations, the significance level of the
%   Rule N test may be slightly higher than the user-specified value). 
%   *** NOTE: To ensure that   p = sig_p,  choose values of p and MC such 
%   that p*MC is an integer. Equivalently, MC must be a multiple of 1/p ***
%
%   metadata: Information to support replication of the analysis.
%
%
% ----- Written By -----
% 
% Jonathan King, 2018, University of Arizona (jonking93@email.arizona.edu)

% Parse the inputs
[plotting, algorithm, sigTest, mcssaArgs, p, convergeTest] = setup(varargin{:});

% Declare the initial structure
s = struct();

% Run an SSA on the time series
[s.singVals, s.singVecs, s.tsm0, s.T, s.C] = simpleSSA(ts, M, algorithm);

% Get the RCs
s.RCs = getRCs( s.singVecs, s.T, algorithm);

% If testing significance
if sigTest
    % Generate the surrogate singular values
    tic
    [s.surrVals] = MC_SSA(ts, M, algorithm, s.singVecs, mcssaArgs{:} );
    toc
    
    % Get the upper level of the 2-tailed significance test
    p2tail = p/2;
    
    % Do a significance test using the surrogate singular values
    [sigThreshold, s.sig_p] = mcSigThreshold( s.surrVals, p2tail);
    s.sigThreshold = permute( sigThreshold, [2,1,3]);
    
    % Get the indices of significant singular values
    s.isSigVal = squeeze( s.singVals >= s.sigThreshold );
        
    % Record Monte Carlo convergence data
    if convergeTest
        [s.MCsigThresh, s.MC_p] = mcSigThreshold(s.surrVals, p2tail, 'convergeTest');
    end
end

% Estimate the periods/frequencies of the singular vectors
s.singPeriod = singVecPeriod(s.singVecs);
s.singFreq = 1 ./ s.singPeriod;

% Assign the metadata
s.metadata = [{'M';'algorithm'}, {M;algorithm}];
if sigTest
    s.metadata(3:5,1:2) = [{'MC';'noise';'p'}, {mcssaArgs{1};mcssaArgs{2};p}];
end

% Make plots if desired
if plotting
    ssaconvergence(s);
    ssasignificance(s);
end

end

%%%%% Helper Functions %%%%%
function[plotting, algorithm, sigTest, mcssaArgs, p, convergeTest] = setup( varargin)

% Parse the inputs
[plotting, algorithm, MC,     noise,    p, sigTest,     parallel, estimateRuntime, showProgress, convergeTest] = parseInputs( varargin,...
{'noplot','algorithm','MC','noiseType','p','noSigTest','parallel','estimateRuntime','showProgress','noConvergeTest'},...
{true,      'BK',    1000,    'red',   0.05, true,    false,        false,            false,         true},...
{ 'b',    {'BK','VG'}, {}, {'red','white'},{}, 'b',       'b',         'b',              'b',           'b'} );

% Get the arguments for MC_SSA
mcssaArgs = {MC, noise};
if parallel
    mcssaArgs = [mcssaArgs, 'parallel'];
end
if estimateRuntime
    mcssaArgs = [mcssaArgs, 'estimateRuntime'];
end
if showProgress
    mcssaArgs = [mcssaArgs, 'showProgress'];
end

end
