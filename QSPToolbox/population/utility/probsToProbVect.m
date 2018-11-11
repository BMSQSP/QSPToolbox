function myProbVect = probsToProbVect(myProbs)
% Reshape a mxn table of axis probilities (m axis, n bins in each axis)
% to a 1 X m*n vector.
% This is needed in order to feed the vector to the optimization
% algorithms
%
% ARGUMENTS
%  myProbs: an array of axis probabilities, preferentially already
%  transformed.
%
    myProbVect = reshape(transpose(myProbs),[],1);
    myProbVect = transpose(myProbVect);
end