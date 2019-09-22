function myVPop = compositeGOF(myVPop)
% Evaluate the composite goodness-of-fit (Fisher's combined test) for a 
% given axes weight solution
%
% ARGUMENTS
%  myVPop:     An instance of a VPop object. The fields should be
%              populated and individual p-values calculated before calling:
%               mnSDTable
%               binTable
%               gofMn
%               gofSD
%               gofBin
%               gofDist
%               gofDist2D
%               gofCor
%
% RETURNS
%  myVpop:     A VPop object with the composite GOF updated.
%

continueFlag = false;
if nargin > 1
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: myVPop.'])
    continueFlag = false;
elseif nargin > 0
    continueFlag = true;
else
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: myVPop.'])
end

if continueFlag
    if sum(ismember(class(myVPop),{'VPop','VPopRECIST'}))<1
        continueFlag = false;
        warning(['Input VPop not recognized in call to ',mfilename,'.'])
    end       
end

if continueFlag
    myMnSDTable = myVPop.mnSDTable;
    myBinTable = myVPop.binTable;
    myDistTable = myVPop.distTable;
    

    if ~isempty(myMnSDTable)
    % We will only include a term in the composite GOF if the corresponding
    % weight is > 0
        indicesMn = find(myMnSDTable{:,'weightMean'} > 0);
        indicesSD = find(myMnSDTable{:,'weightSD'} > 0);
    else
        indicesMn = [];
        indicesSD = [];
    end
    
    if ~isempty(myBinTable)
        indicesBin = find(myBinTable{:,'weight'} > 0);
    else
        indicesBin = [];
    end
    
    if ~isempty(myDistTable)
        indicesDist = find(myDistTable{:,'weight'} > 0);
    else
        indicesDist = [];
    end    
    pVals = [myVPop.gofMn(indicesMn);myVPop.gofSD(indicesSD);myVPop.gofBin(indicesBin);myVPop.gofDist(indicesDist)];
    dof = 2*length(pVals);
    fisherStat = -2*sum(log(pVals));
    myVPop.gof = chi2cdf(fisherStat,dof,'upper');
else
    warning(['Unable to run ',mfilename,'.  Returning input.'])
end    
end