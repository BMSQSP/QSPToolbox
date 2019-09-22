function myDistTable = convertExpDataTo2DDistTable(myVPop,myExpVars1,myExpVars2,mySimTimepoints1,mySimTimepoints2,myInterventions1,myInterventions2)
% This function takes experimental data from a newer 2D table and
% converts it to a 2D dist table format for use in MAPEL
%
% ARGUMENTS:
% myVPop:                 A VPop object with a populated expData field.  A
%                          mapelOptions structure is also OK.
% myExpVars1,myExpVars2:  Two paired cell arrays of length N 
%                          (number of 2D distributions to calibrate to)
%                          with the experimental/data variable names
% mySimTimepoints1,mySimTimepoints2:  Two paired cell arrays of length N 
%                                     (number of 2D distributions to calibrate to)
%                                     with the simulation timepoints
% myInterventions1,myInterventions2:  Two paired cell arrays of length N 
%                                     (number of 2D distributions to calibrate to)
%                                      with the interventions
%
%
% RETURNS
% myDistTable
%

continueFlag = true;
if nargin > 7
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: myVPop (or mapelOptions),myExpVars1,myExpVars2,mySimTimepoints1,mySimTimepoints2,myInterventions1,myInterventions2.'])
    continueFlag = false;
elseif nargin < 7
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: myVPop (or mapelOptions),myExpVars1,myExpVars2,mySimTimepoints1,mySimTimepoints2,myInterventions1,myInterventions2.'])
    continueFlag = false;
else
    continueFlag = true;
end

if continueFlag
    if sum(ismember({'VPop','VPopRECIST','VPopRECISTnoBin','mapelOptions','mapelOptionsRECIST','mapelOptionsRECISTnoBin'},class(myVPop))) < 1
        warning(['Wrong input arguments for ',mfilename,'. Should provide a RECIST type myVPop (or mapelOptions).'])
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
    nRows = length(myExpVars1);
    if isa(myVPop,'VPop') || isa(myVPop,'mapelOptions')
        nDataHeaderCols = 7;
    else
        nDataHeaderCols = 11;
    end
    for rowCounter = 1 : nRows
        myExpVars1Cur = myExpVars1{rowCounter};
        myExpVars2Cur = myExpVars2{rowCounter};
        mySimTimepoints1Cur = mySimTimepoints1(rowCounter);
        mySimTimepoints2Cur = mySimTimepoints2(rowCounter);
        myInterventions1Cur = myInterventions1{rowCounter};
        myInterventions2Cur = myInterventions2{rowCounter};
        if rowCounter == 1
            % The new headers will be a little different since we are
            % 2 measures on each row
            % This is not yest supported as we will have yet to add
            % PatientIDVar for non-RECIST
            if isa(myVPop,'VPop') || isa(myVPop,'mapelOptions')
                tableVariableNames = {'time1','time2','interventionID1','interventionID2','elementID1','elementID2','elementType1','elementType2','expDataID1','expDataID2','expTimeVarID1','expTimeVarID2','expVarID1','expVarID2'};
            else
                tableVariableNames = {'time1','time2','interventionID1','interventionID2','elementID1','elementID2','elementType1','elementType2','expDataID1','expDataID2','expTimeVarID1','expTimeVarID2','expVarID1','expVarID2','PatientIDVar1','PatientIDVar2','TRTVar1','TRTVar2','BRSCOREVar1','BRSCOREVar2','RSCOREVar1','RSCOREVar2'};
            end
            oldTableVariableNames = myVPop.expData.Properties.VariableNames(1:nDataHeaderCols);
            tableVariableNames = [tableVariableNames,{'weight','expN', 'expSample', 'predN', 'predIndices', 'predSample', 'predProbs'}];
            myDistTable = cell2table(cell(0,length(tableVariableNames)));
            myDistTable.Properties.VariableNames = tableVariableNames;
        end
        sourceRow1 = find(ismember(myVPop.expData.('expVarID'),myExpVars1Cur) & ismember(myVPop.expData.('time'),mySimTimepoints1Cur) & ismember(myVPop.expData.('interventionID'),myInterventions1Cur));
        sourceRow2 = find(ismember(myVPop.expData.('expVarID'),myExpVars2Cur) & ismember(myVPop.expData.('time'),mySimTimepoints2Cur) & ismember(myVPop.expData.('interventionID'),myInterventions2Cur));
        curData = myVPop.expData{[sourceRow1,sourceRow2],nDataHeaderCols+1:end};
        curData = curData(:,find(min(~isnan(curData),[],1)));
        %[~,idx] = sort(curData(1,:),'ascend');
        %curData=curData(:,idx);
        [~,nEntries] = size(curData);
        if isa(myVPop,'VPop') || isa(myVPop,'mapelOptions')
            curRow = {myVPop.expData{sourceRow1,'time'},myVPop.expData{sourceRow2,'time'},myVPop.expData{sourceRow1,'interventionID'},myVPop.expData{sourceRow2,'interventionID'},myVPop.expData{sourceRow1,'elementID'},myVPop.expData{sourceRow2,'elementID'},myVPop.expData{sourceRow1,'elementType'},myVPop.expData{sourceRow2,'elementType'},myVPop.expData{sourceRow1,'expDataID'},myVPop.expData{sourceRow2,'expDataID'},myVPop.expData{sourceRow1,'expTimeVarID'},myVPop.expData{sourceRow2,'expTimeVarID'},myVPop.expData{sourceRow1,'expVarID'},myVPop.expData{sourceRow2,'expVarID'},myVPop.expData{sourceRow1,'PatientIDVar'},myVPop.expData{sourceRow2,'PatientIDVar'}};
        else
            curRow = {myVPop.expData{sourceRow1,'time'},myVPop.expData{sourceRow2,'time'},myVPop.expData{sourceRow1,'interventionID'},myVPop.expData{sourceRow2,'interventionID'},myVPop.expData{sourceRow1,'elementID'},myVPop.expData{sourceRow2,'elementID'},myVPop.expData{sourceRow1,'elementType'},myVPop.expData{sourceRow2,'elementType'},myVPop.expData{sourceRow1,'expDataID'},myVPop.expData{sourceRow2,'expDataID'},myVPop.expData{sourceRow1,'expTimeVarID'},myVPop.expData{sourceRow2,'expTimeVarID'},myVPop.expData{sourceRow1,'expVarID'},myVPop.expData{sourceRow2,'expVarID'},myVPop.expData{sourceRow1,'PatientIDVar'},myVPop.expData{sourceRow2,'PatientIDVar'},myVPop.expData{sourceRow1,'TRTVar'},myVPop.expData{sourceRow2,'TRTVar'},myVPop.expData{sourceRow1,'BRSCOREVar'},myVPop.expData{sourceRow2,'BRSCOREVar'},myVPop.expData{sourceRow1,'RSCOREVar'},myVPop.expData{sourceRow2,'RSCOREVar'}};
        end 
        curRow = [curRow,{1, nEntries, {curData}, nan, {nan}, {nan}, {nan}}];
        curRow = cell2table(curRow);
        curRow.Properties.VariableNames = myDistTable.Properties.VariableNames; 
        myDistTable = [myDistTable; curRow];   
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])
end            

end