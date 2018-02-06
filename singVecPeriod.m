function[per] = singVecPeriod( singVec )
%% Approximates the period of a singular vector
%
% Each column is a singular vector

if ~ismatrix(singVec)
    error('singVec must be a matrix');
end
if isrow(singVec)
    singVec = singVec';
end

% Preallocate output
nVecs = size(singVec,2);
per = NaN(nVecs,1);

% For each singular vector
for k = 1:nVecs

    % Get all the local maxima
    [~,localMax] = findpeaks( singVec(:,k) );

    % Get all the local minima
    [~, localMin] = findpeaks( -singVec(:,k) );

    % Sort everything into a single array
    locals = sort([localMax; localMin]);


    % Calculate the half period distances
    nLoc = numel(locals);
    if nLoc > 1
        dist = locals(2:end) - locals(1:end-1);
        per(k) = 2*mean(dist);
    elseif nLoc==1
        per(k) = 2 * length(singVec(:,k));
    end
end

end