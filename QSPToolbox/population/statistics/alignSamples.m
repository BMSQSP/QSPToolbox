function [sample1Ind, sample2Ind, SC] = alignSamples(sample1, sample2)
% This function takes two samples, combines them, and returns
% indices vectors that can be used to map the
% original values in each sample onto the combined vector, prioritizing
% the last observation to facilitate generating empirical CDFs.
%
% ARGUMENTS
%  sample1     observed values for sample1. These are assumed to be sorted
%              in ascending order
%  sample2     observed values for sample2. These are assumed to be sorted
%              in ascending order    
%
% RETURNS
%  sample1Ind  a 1xlength unique sample indices that map elements in sample
%              1 to SC.
%  sample2Ind  a 1xlength unique sample indices that map elements in sample
%              2 to SC.
%  SC          combined sample values
%
	% We need to re-calculate the cdf for sample1 and sample2
	% at all x-values in either
	% sample 1 and sample2 so we can directly calculate
	% the differences in the cdf at each point
    [sample1,Inu1,~]=unique(sample1,'last');
    [sample2,Inu2,~]=unique(sample2,'last');
    SC = unique([sample1,sample2]);
    SC = sort(SC, 'ascend');
    sample1Ind = interp1(sample1,[1:length(sample1)],SC,'previous');
    sample2Ind = interp1(sample2,[1:length(sample2)],SC,'previous');
    % Right fill with last value, leave the left size as nan
    if sum(isnan(sample1Ind)) > 0
        t = 1:numel(sample1Ind);
        sample1Ind = interp1(t(~isnan(sample1Ind)),sample1Ind(~isnan(sample1Ind)),t,'previous','extrap');    
    end
    if sum(isnan(sample2Ind)) > 0
        t = 1:numel(sample2Ind);
        sample2Ind = interp1(t(~isnan(sample2Ind)),sample2Ind(~isnan(sample2Ind)),t,'previous','extrap');    
    end   
    % To update CDFs, we will need to map back to the non-unique but sorted
    % samples
    curInd = find(~isnan(sample1Ind));
    sample1Ind(curInd) = Inu1(sample1Ind(curInd));
    curInd = find(~isnan(sample2Ind));
    sample2Ind(curInd) = Inu2(sample2Ind(curInd));
end