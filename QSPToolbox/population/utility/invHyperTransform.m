function probVect = invHyperTransform(probVectTrans)
% Recover probabilities from hyperspherical coordinate transform
%
% ARGUMENTS
% probVectTrans:      the transformed probability vector
%
% RETURNS
% probVect:           a vector of probabilities (0-1, sum to 1)
%

% We likely want a row vector out of the function
nProbs = length(probVectTrans) + 1;
probVect = ones(1, nProbs);

probVect(1:(nProbs-1)) = probVect(1:(nProbs-1)).*cos(probVectTrans);
for i = 2:nProbs
    probVect(i:nProbs) = probVect(i:nProbs)*sin(probVectTrans(i-1));
end
probVect = probVect.^2;
end