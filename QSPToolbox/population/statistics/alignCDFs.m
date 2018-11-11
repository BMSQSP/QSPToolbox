function [CDF1, CDF2] = alignCDFs(sample1, sample2, W1, W2)
% This function takes two samples and weights, and returns
% adjusted sample and CDF vectors that can be compared element-wise
%
% ARGUMENTS
%  sample1     observed values for sample1. These are assumed to be sorted
%              in ascending order
%  sample2     observed values for sample2. These are assumed to be sorted
%              in ascending order    
%  W1          a 1xlength(sample1) vector of weights for sample 1
%  W2          a 1xlength(sample2) vector of weights for sample 2
%
% RETURNS
%  CDF1        a 1xlength unique sample values of increasing cumulative
%              probabilities in sample 1 corresponding to the unique 
%              combined sample values
%  CDF2        a 1xlength unique sample values of increasing cumulative
%              probabilities in sample 1 corresponding to the unique 
%              combined sample values
% First we need to interpolate the two samples so we can directly
% take the supremum function.
    W1 = W1 / sum(W1);
    W2 = W2 / sum(W2);
    W1 = cumsum(W1);
    W2 = cumsum(W2);
	% We need to re-calculate the cdf for sample1 and sample2
	% at all x-values in either
	% sample 1 and sample2 so we can directly calculate
	% the differences in the cdf at each point
    [sample1ind, sample2ind, SC] = alignSamples(sample1, sample2);
    CDF1 = zeros(1, length(SC));
    CDF2 = zeros(1, length(SC));
    CDF1Indices = find(~isnan(sample1ind));
    CDF2Indices = find(~isnan(sample2ind));
    CDF1(CDF1Indices) = W1(sample1ind(CDF1Indices));
    CDF2(CDF2Indices) = W2(sample2ind(CDF2Indices));
    
end