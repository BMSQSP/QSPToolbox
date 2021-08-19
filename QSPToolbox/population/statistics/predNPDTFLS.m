function myVPop = predNPDTFLS(myVPop)
% This function takes a VPop object and updates with 
% predictions for number of PD patients.
%
% ARGUMENTS
%  myVPop:                  (Required) a mapelOptionsRECIST or VPopRECIST
%                                     instance.
%
% Returns
%  myVPop:                  VPopRECIST with updated brTableRECIST.  This 
%                            gives predNPD21LS values.
%                            
%


continueFlag = true;
if nargin > 1
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: VPopRECIST object.'])
    continueFlag = false;
elseif nargin < 1
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: VPopRECIST object.'])
end

if continueFlag
    if ~ismember(class(myVPop),{'VPopRECIST'})
        continueFlag = false;
        warning(['Input myVPop not recognized in call to ',mfilename,'.  Requires: VPopRECIST object.'])
    end       
end


if continueFlag
    lmResults = regressPD2ScaleFactor(myVPop);
    % Update the values for the 
    % First find the blocks by the PD2 measures
    testRows = find(sum(isnan(myVPop.brTableRECIST{:,{'predN','predCR','predPR','predSD','predPD'}}),2)==0);
    % Just use the first time point for the PD2
    testData = myVPop.brTableRECIST(:,{'interventionID'});
    [C,IA,IC] = unique(testData,'rows','stable');    
    predRows = intersect(testRows, IA);
    predVals = myVPop.brTableRECIST{predRows,{'predSD','predPD'}};
    predN = myVPop.brTableRECIST{predRows,{'predN'}};
    res = predVals*lmResults.Coefficients{:,'Estimate'}.*predN;
    % Maybe there is a faster assignment, but scan through the
    % the rows for now.
    for rowCounter = 1:length(testRows)
        curRow = testRows(rowCounter);
        if ismember(curRow,testRows) > 0
            sourcePred = max(find(IA<=curRow));
            myVPop.brTableRECIST{rowCounter,'predNPD21LS'}=res(sourcePred);
        end
    end
else
    lmResults = {};
    warning(['Unable to run ',mfilename,'.'])
end
end


