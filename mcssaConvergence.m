function[MCsigThresh, MC_p, MCisSig] = mcssaConvergence( surrVals, p, singVals)

% Run a convergence test on the Monte Carlo iterations
[MCsigThresh, ~, MC_p] = mcthreshold( surrVals, p, '2tail', 'converge');

% Preallocate a logical matrix for whether a given singular value is
% significant at a particular iteration.
nMC = size(MCsigThresh,1);
MCisSig = false( size(MCsigThresh) );

% Get the development of significant RCs
for k = 1:nMC
    MCisSig(k,:) = singVals > MCsigThresh(k,:);
end

end
    
 