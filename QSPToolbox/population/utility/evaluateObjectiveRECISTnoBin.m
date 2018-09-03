function myObjective = evaluateObjectiveRECISTnoBin(myVPop, myPWVectTrans)
% This function evaluates an objective function, which is similar to
% the composite goodness-of-fit but with a few critical differences, and is
% generally called only during the optimization / VPop fit to data.
%
% ARGUMENTS
%  myVPop:              a VPopRECISTnoBin object instance, with properties already
%                       assigned:
%                        binProbs:  Should at least be initialized
%                        spreadOut
%                        minEffN
%                        mnSDTable: should at least be populated with
%                                   experiment data
%                        binTable:  Also need experiment data
%  myPWVectTrans:     a 1 x (nVP-1) vector of transformed PW
%                       probabilities
%
% RETURNS
%  myObjective:         the objective function values
%

% Next update the individual GOF statistics
myVPop.pws=invHyperTransform(myPWVectTrans);
myVPop = evaluateGOF(myVPop);

myObjective = 0;
effN = 1/(sum((myVPop.pws).^2));
if myVPop.spreadOut > 0
    nPWs = length(myVPop.pws);
    notUnif = myVPop.spreadOut*(nPWs/effN-1)/(nPWs-1);
    myObjective = myObjective + notUnif;
end

if myVPop.minEffN > 0
    myObjective = myObjective + (effN < myVPop.minEffN) * 1E6 * (myVPop.minEffN - effN);
end

myMnSDTable = myVPop.mnSDTable;
myBinTable = myVPop.binTable;
myDistTable = myVPop.distTable;
myBRTable = myVPop.brTableRECIST;
myRTable = myVPop.rTableRECIST;

% We use the same epsilon from the MAPEL paper
%epsilon = myVPop.epsilon;
epsilon = 1E-16;
gofMean = myVPop.gofMn;
gofSD = myVPop.gofSD;
gofBin = myVPop.gofBin;
gofDist = myVPop.gofDist;
gofBR = myVPop.gofBR;
gofR = myVPop.gofR;

% May want to add checks here for the Mean/SD/Bin evaluations and whether
% there are entries before adding the respective terms
if ~isempty(myMnSDTable)
    myObjective = myObjective - 2*(sum(myMnSDTable{:,'weightMean'} .* log10(gofMean+epsilon)) + sum(myMnSDTable{:,'weightSD'} .* log10(gofSD+epsilon)));
end
if ~isempty(myBinTable)
    myObjective = myObjective - 2*(sum(myBinTable{:,'weight'} .* log10(gofBin+epsilon)));
end
if ~isempty(myDistTable)
    myObjective = myObjective - 2*(sum(myDistTable{:,'weight'} .* log10(gofDist+epsilon)));
end
if ~isempty(myBRTable)
    myObjective = myObjective - 2*(sum(myBRTable{:,'weight'} .* log10(gofBR+epsilon)));
end
if ~isempty(myRTable)
    myObjective = myObjective - 2*(sum(myRTable{:,'weight'} .* log10(gofR+epsilon)));
end
end