function pVal = contingency24(observations1, observations2, exactFlag)
% This function returns the p-value for a two-sided test of whether  
% the null hypothesis that two samples with four categories come from the
% same distribution.
% For small sample sizes, i.e. if any of the observations are < 5,
% an exact test is applied.  For sufficient sample sizes, the chisq  
% approximation is applied.
%
% ARGUMENTS
%  observations1: a 1x4 vector of observations in each bin
%  observations2: a 1x4 vector of observations in each bin
%  exactFlag:     boolean, whether to check if an exact test should be
%                 performed.  If false, the chi-square approximation is
%                 imposed.
%
% RETURNS
%  pVal
%
pVal = nan;
continueFlag = true;

if ((max(size(observations1))~=4) || (min(size(observations1))~=1) || (max(size(observations2))~=4) || (min(size(observations2)~=1)))
    continueFlag = false;
    warning(['Unable to continue in ',mfilename,', two 1x4 vectors of categorical observation counts should be supplied.'])
end

if continueFlag
    no1t = sum(observations1);
    no2t = sum(observations2);
    ep1 = (observations1(1) + observations2(1)) / (no1t + no2t);
    ep2 = (observations1(2) + observations2(2)) / (no1t + no2t);
    ep3 = (observations1(3) + observations2(3)) / (no1t + no2t);
    ep4 = (observations1(4) + observations2(4)) / (no1t + no2t);
    epAll = [ep1, ep2, ep3, ep4];
    % We use the rule of thumb:
    % mentioned in several places, including:
    % McHugh, ML, Biochem Med (Zagreb). 2013 Jun; 23(2): 143–149.
    % No more than 20% of the expected counts are less than 5 and 
    % all individual expected counts are 1 for the chisq approximation
    % to hold
    minFail = min(min([no1t, no2t]) * epAll) < 1;
    fiveFail = sum(sum(([no1t; no2t] * epAll) < 5)) > (0.2*8);
    zeroFail = min(epAll) <= 0;
    
   if ((minFail || fiveFail) & (exactFlag || zeroFail))
       chisqFlag = false;
   else
        chisqFlag = true;
   end
end

if continueFlag
    if chisqFlag
        expected = [ep1*no1t, ep2*no1t, ep3*no1t, ep4*no1t, ep1*no2t, ep2*no2t, ep3*no2t, ep4*no2t];
        % We already apply our own check on the minimum number of observations
        % per bin.  Don't use the default 'emin' of 5.
        [h,pVal,stats] = chi2gof([1 2 3 4 5 6 7 8],'freq',[observations1, observations2],'expected',expected,'ctrs',[1 2 3 4 5 6 7 8],'nparams',4,'emin',0);
    else
      
        x = [observations1; observations2];
  
        % We use the 3rd party function for this calculation
        % There were some issues that we have debugged.
        [T,pVal] = evalc('myfisher24(x);');
    end
else
    warning(['Unable to continue in ',mfilename,', two 1x4 vectors of categorical observation counts should be supplied.'])
end
end