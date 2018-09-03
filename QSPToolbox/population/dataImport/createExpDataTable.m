function myExpDataTable = createExpDataTable(myWorksheet, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID, previousExpDataTable)
% Convert experimental data attached to a worksheet
% into a format that can
% can be used in the population weighting algorithms.
%
% ARGUMENTS
% myWorksheet:           a worksheet containing the data
% myInterventionID:      the ID of the worksheet intervention to map the 
% elementID:             data to the variable ID in the simulation results
% elementType:           type of variable to map to: parameter, 
%                        compartment, species
% expDataID:             dataset ID for the experimental data
% expTimeVarID:          name of the time variable in the experimental 
%                        dataset to map onto the simulation time.  Note if 
%                        the simulated intervention uses a delay or time 
%                        shift, there may need to be a shift in the 
%                        experimental data here to align the two.  Time 
%                        units should be consistent.
% expVarID:              Column ID for the experimental data to map
% previousExpDataTable:  a previous output of the function, if provided the 
%                        additional results will be appended to the end. 
%
% RETURNS
% myExpDataTable:       A MATLAB table with the data and mappings.  This will
%                   have headings:
%                   'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'expVal1', 'expVal2', ... , 'expValN'
%
continueFlag = true;
myExpDataTable = table();
if nargin > 8
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: myWorksheet, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID; optionally: previousExpDataTable.'])
    continueFlag = false;
elseif nargin > 7
    continueFlag = true;
elseif nargin > 6
    previousExpDataTable = table();
else
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: myWorksheet, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID; optionally: previousExpDataTable.'])
    continueFlag = false;
end

if continueFlag
    worksheetInterventionIDs = getInterventionIDs(myWorksheet);
    if sum(ismember(worksheetInterventionIDs,myInterventionID)) < 1
        continueFlag = false;
        warning(['Missing specified InterventionID in worksheet supplied to ',mfilename,'.'])
    end
    if sum(ismember(fields(myWorksheet),'compiled')) < 1
        continueFlag = false;
        warning(['Worksheet supplied to ',mfilename,' should have a compiled model.'])
    else
        elementType = lower(elementType);
        myWorksheetElements = myWorksheet.compiled.elements;
        if sum(ismember(myWorksheetElements(:,2),elementType)) < 1
            continueFlag = false;
            warning(['The elementType supplied to ',mfilename,' should be species, compartment, or parameter.'])
        else
            if sum(ismember(myWorksheetElements(:,1),elementID)) < 1
                continueFlag = false;
                warning(['The elementID supplied to ',mfilename,' was not identified in the compiled model.'])
            end
        end
    end
    worksheetExpDataIDs = getExpDataIDs(myWorksheet);
    if sum(ismember(worksheetExpDataIDs,expDataID)) < 1
        continueFlag = false;
        warning(['The expDataID supplied to ',mfilename,' was not identified in the supplied worksheet.']);       
    else
        myExpTable = getExpData(expDataID, myWorksheet);
        if sum(ismember(myExpTable.Properties.VariableNames, expTimeVarID)) < 1
            continueFlag = false;
            warning(['The expTimeVarID supplied to ',mfilename,' was not identified in the specified dataset.'])
        end
        if sum(ismember(myExpTable.Properties.VariableNames, expVarID)) < 1
            continueFlag = false;
            warning(['The expVarID supplied to ',mfilename,' was not identified in the specified dataset.'])
        end
            
    end
end
  
if continueFlag
    tableVariableNamesFixed = {'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID'};
    [nPreviousRows, nPreviousColumns] = size(previousExpDataTable);
    if nPreviousRows >= 1
        previousVariableNames = previousExpDataTable.Properties.VariableNames;
        diffNames = setdiff(previousVariableNames,tableVariableNamesFixed);        
        if sum(ismember(tableVariableNamesFixed, previousVariableNames)) ~= length(tableVariableNamesFixed)
            continueFlag = false;
            warning(['When providing a previousExpDataTable to ',mfilename,', the VariableNames must include: ',strjoin(tableVariableNamesFixed,', '),'.'])
        elseif sum(strncmpi('expVal',diffNames,6)) ~= length(diffNames)
            continueFlag = false;
            warning(['When providing a previousExpDataTable to ',mfilename,', the VariableNames must include only: ',strjoin(tableVariableNamesFixed,', '),', and "expValXX".'])       
        else
            existingN = length(diffNames);
            myExpDataTable = previousExpDataTable;
        end
    else
        existingN = 0;
        myExpDataTable = cell2table(cell(0,length(tableVariableNamesFixed)));
        myExpDataTable.Properties.VariableNames = tableVariableNamesFixed;
    end
end

if continueFlag
    expDataTime = myExpTable.(expTimeVarID);
    expDataVar = myExpTable.(expVarID);
    the_indices = find((~isnan(expDataTime)) & (~isnan(expDataVar)));
    expDataTime = expDataTime(the_indices);
    expDataVar = expDataVar(the_indices);
    uniqueExpDataTime = unique(expDataTime);
    for timePointCounter = 1 : length(uniqueExpDataTime)
        curTime = uniqueExpDataTime(timePointCounter);
        originalRowIndices = find(expDataTime == curTime);
        allValues = (expDataVar(originalRowIndices));
        curN = length(originalRowIndices);
        if curN > existingN
            newNames = cell(1,0);
            for newNCounter = (existingN + 1) : curN
                newNames = [newNames,['expVal',num2str(newNCounter)]];
            end
            [curNRows, curNCols] = size(myExpDataTable);
            appendTable = nan(curNRows,curN-existingN);
            appendTable = array2table(appendTable,'VariableNames',newNames);
            myExpDataTable = [myExpDataTable, appendTable];
            existingN = curN;
            curRow = [{curTime, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID},num2cell(transpose(allValues))];
        elseif curN < existingN
            paddingValues = nan(1,existingN-curN);
            curRow = [{curTime, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID},num2cell(transpose(allValues)),num2cell(paddingValues)];
        else
            curRow = [{curTime, myInterventionID, elementID, elementType, expDataID, expTimeVarID, expVarID},num2cell(transpose(allValues))];
        end
        curRow = cell2table(curRow);
        curRow.Properties.VariableNames = myExpDataTable.Properties.VariableNames; 
        myExpDataTable = [myExpDataTable; curRow];
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])
end

end