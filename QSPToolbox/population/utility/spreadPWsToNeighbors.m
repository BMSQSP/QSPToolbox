function myNewPWs = spreadPWsToNeighbors(myPWs,myCoefficients,nNew)
% This is a "utility function" that should be called directly.
% This function will generate new PW guesses spread
% around an initial guess, with weights given
% preferentially to nearer neighbors.
%
% ARGUMENTS
%  myPWs:          			   a 1 x nVPs vector of prevalence weights
%  myCoefficients			   an nAxis x nVP matrix of coefficients
%  nNew                        number of new PW solutions to add
%
% RETURNS
%  myNewPWs:          		   original PWs appended with additional
%                              (nNew + 1) x nVP matrix
%

[nDim,nParents] = size(myCoefficients);
vpDistance = pdist2(myCoefficients',myCoefficients');
spreadMuMax = min(4,nParents/2);

% We will generate new guesses that follow exponential
% distributions around the initial point
myNewPWs = zeros(nNew,nParents);
% loop over VPs in the outer loop so we only have to sort distances
% for each VP once
for parentCounter = 1 : nParents
	[~, sortIndex] = sort(vpDistance(parentCounter,:),'ascend');
	lambda = rand([nNew,1])*spreadMuMax;
	for distCounter = 1 : nNew
		spreadVals = exppdf([0:1:(nParents-1)],lambda(distCounter));
		spreadVals = spreadVals./sum(spreadVals);
		myNewPWs(distCounter, sortIndex) = myPWs(parentCounter)*spreadVals+myNewPWs(distCounter, sortIndex);
	end
end
myNewPWs = [myPWs;myNewPWs];
end