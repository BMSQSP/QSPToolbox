function objectiveValue = evaluatePointObjectiveWtRangeCE(uniqueSimTime, simY, expTime, expY)
% Objective function for point response type elements
% WtRange - values being included in the objective are weighted by how far 
%           they fall outside the observed range in the data
% CE - Objective values are also normalized by an error term that is 
%      constant across time points, the average standard deviation in the 
%      data
%
% ARGUMENTS
%  uniqueSimTime: filtered to match unique values in expTime
%  simY: interpolated simulation values at each uniqueSimTime
%  expTime: experiment times for samples
%  expY: experimental Y values
%
% RETURNS
%  objectiveValue: a single value
%
nExpPoints = length(expY);
objectiveValue = 0;
nTpoints = length(uniqueSimTime);
rangeNorm = 0;
for timeCounter = 1 : nTpoints
    curTime = uniqueSimTime(timeCounter);
    curY = simY(timeCounter);
    expIndices = find(expTime == curTime);
    nCurPoints = length(expIndices);
    maxval = max(expY(expIndices));
    minval = min(expY(expIndices));
    % We could normalize out of bounds
    % simulation points by stdev, range, etc...
    % just to provide a gentle push to bring these in
    % spec
    % rangeNorm = (maxval - minval);
    rangeNorm = rangeNorm + std(expY(expIndices));
    flagTop = curY > maxval;
    flagBottom = curY < minval;
    % This was set to 0.99 prior to 12/26/15, then tried adjusting
    % to see if I could get more agreeable results from the optimizers.
    % inRangeWeight = 0.0 seems to give better results for many
    % scenarios
    inRangeWeight = 0.0;
    curVal = inRangeWeight*((flagTop)+(flagBottom))+(1-inRangeWeight)*(flagBottom*(minval-curY)+flagTop*(curY-maxval));
    %curVal = sqrt((sum((expY(expIndices) - curY).^2))/nCurPoints)  / abs(mean(expY(expIndices)));
    objectiveValue = objectiveValue + curVal*1/nTpoints;
end
rangeNorm = rangeNorm / nTpoints;
objectiveValue = objectiveValue / rangeNorm;
end