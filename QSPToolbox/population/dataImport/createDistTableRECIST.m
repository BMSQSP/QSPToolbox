function myDistTableRECIST = createDistTableRECIST(myWorksheet,myExpDataIDs,PatientIDVar,TRTVar,BRSCOREVar,RSCOREVar,myVars,mySimVars,mySimVarTypes,timeVar,startTime,myInterventionIDs)
% This function takes experimental data embedded in a worksheet and
% converts it to a distribution table, bypassing reading from a mapelOptions
% "experimental data table" and directly taking data from the worksheet,
% accounting for whether patients should be excluded because of prior
% CR, PD, or they are clearly of treatment.
%
% Note the cell arrays below are based on calibration variable x interventionID
% combinations.  That is, if multiple time points are available they will automatically be extracted.
% You can go back and remove rows you do not want or set weights to 0.
%
% ARGUMENTS:
%  myWorksheet:       A worksheet with the experimental data attached
%  myExpDataIDs:      A cell array, of length N intervention/var combinations
%                      with the experimental dataset IDs
%  PatientIDVar:      ID for the variable with the patient IDs
%  TRTVar:            ID for the binary variable indicating whether 
%                      a patient is on treatment
%  BRSCOREVar:        ID for the best overall response to date variable.
%  RSCOREVar:         ID for the current response variable.
%  myVars:            A cell array, of length N intervention/var combinations
%                      variables in the dataset to get mean/SD summary data for.
%  mySimVars:         Simulation variable names, of length N intervention/var
%                      combinations, map to myVars.
%  mySimVarTypes:     A cell array, of length N intervention/var combinations
%                      with values to indicate the type of simulation
%                      variables, i.e.
%                      'parameter', 'compartment', 'species'
%  timeVar:           ID for the variable with time in the experimental dataset
%                     NOTE: it is assumed start of therapy is t=startTime.
%  startTime		  Time to start therapy.
%  myInterventionIDs: A cell array, of length N intervention/var combinations
%                      with all of the intervention names from the simulations
%
% RETURNS
%  myDistTableRECIST

% TODO: Add proofing of inputs
continueFlag = true;

tableVariableNamesFixed = {'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'PatientIDVar', 'TRTVar', 'BRSCOREVar', 'RSCOREVar'};
tableVariableNames = [tableVariableNamesFixed,{'weight','expN', 'expSample', 'predN', 'predIndices','predSample', 'predProbs','expCombinedIndices','simCombinedIndices','combinedPoints'}];
myDistTableRECIST = cell2table(cell(0,length(tableVariableNames)));
myDistTableRECIST.Properties.VariableNames = tableVariableNames;

% Proofing of the input
% variables starts here
if continueFlag
    
    allExpDataIDs = getExpDataIDs(myWorksheet);
    allInterventionIDs = getInterventionIDs(myWorksheet);
    
    for varCounter = 1 : length(myVars)
    
        curVar = myVars{varCounter};
        myExpDataID = myExpDataIDs{varCounter};
        interventionID = myInterventionIDs{varCounter};
  
        curIndex = find(ismember(allExpDataIDs,myExpDataID));
        
        if isempty(curIndex)
            warning(['Not able to find experimental dataset ',myExpDataID,' in ',mfilename,'.'])
            continueFlag = false;
        else        
            curData = myWorksheet.expData{curIndex}.data;
            % First remove off treatment rows
            curData = curData(find(curData{:,TRTVar}==1),:);
            [nRows,nCols] = size(curData);
            if nRows < 1
                warning(['Not able to find on treatment data for experimental dataset ',myExpDataID,' in ',mfilename,'.'])
                continueFlag = false;
            else
                if sum(ismember(curData.Properties.VariableNames,myVars{varCounter})) < 1
                    warning(['Not able to find experimental variable ',myVars{varCounter},' in dataset ',myExpDataID,' in ',mfilename,'.'])
                    continueFlag = false;
                end
                if sum(ismember(allInterventionIDs,interventionID)) < 1
                    warning(['Not able to find intervention ID ',interventionID,' associated with variable ',myVars{varCounter},' in dataset ',myExpDataID,' in ',mfilename,'.'])
                    continueFlag = false;                
                end
            end
        end
    end
end
     

if continueFlag
    for varCounter = 1 : length(myVars)
        curVar = myVars{varCounter};
        myExpDataID = myExpDataIDs{varCounter};
        interventionID = myInterventionIDs{varCounter};
        curIndex = find(ismember(allExpDataIDs,myExpDataID));
        curData = myWorksheet.expData{curIndex}.data;        
        
        % First remove off treatment rows
        curData = curData(find(curData{:,TRTVar}==1),:);
        % Need to make all table vars into string for this to work
        new_variable = arrayfun(@(value) cast(value{1}, 'char'), curData.(BRSCOREVar), 'uniform', 0);
        curData.(BRSCOREVar) = new_variable;
        new_variable = arrayfun(@(value) cast(value{1}, 'char'), curData.(RSCOREVar), 'uniform', 0);
        curData.(RSCOREVar) = new_variable;    
        % Now also filter once a patient hits CR, PD
        curPatientIDs = unique(curData{:,PatientIDVar},'stable');
        nPatients = length(curPatientIDs);
        allRows = nan(0,1);
        for patientCounter = 1 : nPatients
            curRows = find(ismember(curData{:,PatientIDVar},curPatientIDs(patientCounter)));
            % First filter out PD2 and subsequent time points
			% We will exclude all data with PD2 for now
			lastRows = find(ismember(curData{curRows,BRSCOREVar},{'PD2'}) | ismember(curData{curRows,RSCOREVar},{'PD2'}));
			if ~isempty(lastRows)
				lastRows = lastRows(1)-1;
                if lastRows > 0
                    curRows = curRows(1:lastRows);
                else
                    curRows = [];
                end
			end               
            % We will try to find where t == 0.  In some cases, such as
            % when we calibrate specifically to post-treatment data,
            % we may have selected out the baseline data to calibrate
            % separately.  In this case, we find the first positive time.
            if length(curRows) > 0            
                firstRow = find(curData{curRows,timeVar}==0);
                if isempty(firstRow)
                    firstRow = find(curData{curRows,timeVar} > 0);
                    firstRow = firstRow(1);
                end            
                lastRows = find(ismember(curData{curRows,BRSCOREVar},{'PD','CR'}) | ismember(curData{curRows,RSCOREVar},{'PD','CR'}));
                otherMeasureRows = find(ismember(curData{curRows,RSCOREVar},{'PR','SD'}));
                if (length(lastRows) > 0)
                    % We will include up to the first CR/PD measure
                    lastRows = lastRows(1);
                elseif length(otherMeasureRows) > 0
                    % Otherwise, we take the last time point where we have a valid
                    % response measure
                    lastRows = otherMeasureRows(length(otherMeasureRows));
                else % if (length(lastRows) < 1) % This should be true...
                    lastRows = firstRow;
                end
                curRows = curRows(firstRow:lastRows);
                allRows = [allRows; curRows];
            end
        end
        selectData = curData(allRows,:);  
    
        curRows = find(~isnan(selectData{:,curVar}));        
        if length(curRows) >= 2
            % Now we need to step through and evaluate for the times
            curTimes = unique(selectData{curRows,timeVar});
			
            for curTimeCounter = 1: length(curTimes)
                curDataRows = find(selectData{curRows,timeVar} == curTimes(curTimeCounter));
                if length(curDataRows) >= 2
                    curData = (selectData{curRows(curDataRows),curVar})';
                    % We only keep add a row if the largest value is bigger
                    % than the smallest.  Sometimes when normailizing to
                    % baseline variables with zero variability are created
                    % at T=0 and it makes sense to filter these.
                    if max(curData) > min(curData)					
                        newRow = {curTimes(curTimeCounter), interventionID, mySimVars{varCounter}, mySimVarTypes{varCounter}, myExpDataID, timeVar, curVar, PatientIDVar, TRTVar, BRSCOREVar, RSCOREVar};
                        newRow = [newRow,{1, length(curData), {sort(curData,'ascend')}, nan, {nan}, {nan}, {nan},{nan},{nan},{nan}}];
                        newRow = cell2table(newRow);
                        %newRow = table(curTimes(curTimeCounter), interventionID, mySimVars{varCounter}, mySimVarTypes{varCounter}, myExpDataID, timeVar, curVar,PatientIDVar,TRTVar,BRSCOREVar,RSCOREVar1, length(curData), {sort(curData,'ascend')}, nan, {nan}, {nan}, {nan},{nan},{nan},{nan})
                        newRow.Properties.VariableNames = myDistTableRECIST.Properties.VariableNames; 
                        myDistTableRECIST = [myDistTableRECIST; newRow];
                    end
                end
            end
        end
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])    
end    