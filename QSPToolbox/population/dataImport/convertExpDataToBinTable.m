function myBinTable = convertExpDataToMnSDTable(myVPop)
% This function takes experimental data and
% converts it to a bin table format for use in MAPEL
%
% ARGUMENTS:
% myVPop:           A VPop object with a populated expData field.  A
%                   mapelOptions structure is also OK.
%
% RETURNS
% myBinTable
%

continueFlag = true;
if nargin > 1
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: myVPop.'])
    continueFlag = false;
elseif nargin > 0
    continueFlag = true;
else
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: myVPop.'])
    continueFlag = false;
end

if continueFlag
    if sum(ismember({'VPop','VPopRECIST','VPopRECISTnoBin','mapelOptions','mapelOptionsRECIST','mapelOptionsRECISTnoBin'},class(myVPop))) < 1
        warning(['Wrong input arguments for ',mfilename,'. Should provide: myVPop (or mapelOptions).'])
        continueFlag = false;
    end
end
        
if continueFlag
    if sum(ismember({'table'},class(myVPop.expData))) < 1
        warning(['Wrong input arguments for ',mfilename,'. Should provide: myVPop (or mapelOptions) with a populated expData property.'])
        continueFlag = false;
    end
end

if continueFlag
    [nRows, ~] = size(myVPop.expData);
    if isa(myVPop,'VPop') || isa(myVPop,'mapelOptions')
        nDataHeaderCols = 8;
    else
        nDataHeaderCols = 11;
    end
    uniqueExpDataTime = unique(myVPop.expData{:,'time'});
    for rowCounter = 1 : nRows
        if rowCounter == 1
            tableVariableNames = myVPop.expData.Properties.VariableNames(1:nDataHeaderCols);
            tableVariableNames = [tableVariableNames,{'weight','binEdge1','binEdge2','binEdge3','expN','expBin1','expBin2','expBin3','expBin4','predN','predBin1','predBin2','predBin3','predBin4'}];
            myBinTable = cell2table(cell(0,length(tableVariableNames)));
            myBinTable.Properties.VariableNames = tableVariableNames;
        end
        curData = myVPop.expData{rowCounter,nDataHeaderCols+1:end};
        curData = curData(~isnan(curData));
        curRow = table2cell(myVPop.expData(rowCounter,1:nDataHeaderCols));
        myElementID = myVPop.expData{rowCounter,'elementID'};
        allCurVarRows = find(ismember(table2cell(myVPop.expData(:, 'elementID')),myElementID));
        allCurVarValues = myVPop.expData{allCurVarRows, nDataHeaderCols+1:end};
        allCurVarValues = allCurVarValues(~isnan(allCurVarValues));
        myBinEdgeValues = [(median(allCurVarValues) - min(allCurVarValues))/2+min(allCurVarValues),median(allCurVarValues),(max(allCurVarValues) - median(allCurVarValues))/2+median(allCurVarValues)];
        curProbs = wtdBinProb(curData, ones(1, length(curData))/length(curData), myBinEdgeValues);
        expN = length(curData);
        curRow = [curRow,{1,myBinEdgeValues(1), myBinEdgeValues(2), myBinEdgeValues(3), expN, curProbs(1), curProbs(2), curProbs(3), curProbs(4), nan, nan, nan, nan, nan}];        
        curRow = cell2table(curRow);
        curRow.Properties.VariableNames = myBinTable.Properties.VariableNames; 
        myBinTable = [myBinTable; curRow];   
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])
end            

end