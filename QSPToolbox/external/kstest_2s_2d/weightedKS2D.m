function myPVal = weightedKS2D(quadrantCounts1cross, quadrantCounts1self, quadrantCounts2cross, quadrantCounts2self, W1, W2, N1, N2)
% This function performs a two-sided 2D KS test and returns a p-value
% It is assumed quadrant counts are performed prior to running this test
%
% References: 
%  Peacock, JA. Two-dimensional goodness-of-fit testing
%  in astronomy, Mon Not R astr Soc 202 (1983)
%  615-627.
%
% This code modified to account for prevalence weighting
% in accordance with license from a MATLAB implementation,
% original written by:
% Dylan Muir (From kstest_2s_2d by Qiuyan Peng @ ECE/HKUST) Date:
% 13th October, 2012
% Modification by Brian Schmidt, August 2018
%
%
% ARGUMENTS
%  quadrantCounts1cross   quadrant counts with each value in sample 1 as reference against sample 2
%  quadrantCounts1self
%  quadrantCounts2cross   quadrant counts with each value in sample 2 as reference against sample 1
%  quadrantCounts2self
%  W1                     a 1xlength(sample1) vector of weights for sample 1
%  W2                     a 1xlength(sample2) vector of weights for sample 2
%  N1                     number of observations underlying sample1
%  N2                     number of observations underlying sample2
%
% RETURNS
%  myPVal                 the p-value resulting from the two-tailed test

nSamples1 = length(W1);
b11 = nSamples1;
b12 = nSamples1*2;
b13 = nSamples1*3;
nSamples2 = length(W2);
b21 = nSamples2;
b22 = nSamples2*2;
b23 = nSamples2*3;
KSstatistic = ones(nSamples1+nSamples2,1);
WQuadC = zeros(nSamples1+nSamples2,4);
WQuadS = zeros(nSamples1+nSamples2,4);

for j = 1:(nSamples1+nSamples2)

   if (j<=nSamples1)
	  WQuadC(j,:) = [sum(quadrantCounts1cross(j,1:b21).*W2),sum(quadrantCounts1cross(j,b21+1:b22).*W2),...
	  sum(quadrantCounts1cross(j,b22+1:b23).*W2),sum(quadrantCounts1cross(j,b23+1:end).*W2)];
	  WQuadS(j,:) = [sum(quadrantCounts1self(j,1:b11).*W1),sum(quadrantCounts1self(j,b11+1:b12).*W1),...
	  sum(quadrantCounts1self(j,b12+1:b13).*W1),sum(quadrantCounts1self(j,b13+1:end).*W1)];
   else
	  j2 = j - nSamples1;
	  WQuadC(j,:) = [sum(quadrantCounts2cross(j2,1:b11).*W1),sum(quadrantCounts2cross(j2,b11+1:b12).*W1),...
	  sum(quadrantCounts2cross(j2,b12+1:b13).*W1),sum(quadrantCounts2cross(j2,b13+1:end).*W1)];
	  WQuadS(j,:) = [sum(quadrantCounts2self(j2,1:b21).*W2),sum(quadrantCounts2self(j2,b21+1:b22).*W2),...
	  sum(quadrantCounts2self(j2,b22+1:b23).*W2),sum(quadrantCounts2self(j2,b23+1:end).*W2)];	  
	  
   end
   KSstatistic(j) = max(abs(WQuadC(j,:) - WQuadS(j,:)));
end 

% Peacock's method for P-value
Neff = N1 * N2 /(N1 + N2);
KSstatistic = max(KSstatistic);
Zn = sqrt(Neff) * KSstatistic;
Zinf = Zn / (1 - 0.53 * Neff^(-0.9));
myPVal = 2 * exp(-2 * (Zinf - 0.5).^2);

% This approximation is not as valid for large
% p-values
if (myPVal > 1)
   myPVal = 1;
end

end