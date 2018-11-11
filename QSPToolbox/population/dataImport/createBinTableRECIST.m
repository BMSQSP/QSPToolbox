function myBinTableRECIST = createBinTableRECIST(myWorksheet,myExpDataIDs,PatientIDVar,TRTVar,BRSCOREVar,RSCOREVar,myVars,myBinEdges,mySimVars,mySimVarTypes,timeVar,startTime,myInterventionIDs)
% This function takes experimental data embedded in a worksheet and
% converts it to a bin table, bypassing the "experimental data table"
% directly accounting for whether patients should be excluded because of
% CR, PD, or they are clearly off treatment.
%
% ARGUMENTS:
%  myWorksheet:       a workshet with the experimental data
%  myExpDataIDs:      a cell array, 1 x n time points and variables
%                      with the experimental data IDs
%  PatientIDVar:      variable with the patient IDs
%  TRTVar:            a binary variable indicating whether a patient is on
%                     treatment
%  BRSCOREVar:        best overall response to date.
%  RSCOREVar:         current response.
%  myVars:            a cell array, 1 x n time points and variables
%                      variables to get mean/SD data for.
%  myBinEdges:    a cell array of length length(myVars),
%                 with each entry containing a cutoff of 3.
%                 if you want to set one atuomatically from data,
%                 just leave the corresponding entry blank. i.e.
%                 {[edge1,edge2,edge3], '', ...}
%  mySimVars:         simulation variable names, map to myVars
%  mySimVarTypes:     a cell array, 1 x n time points and variables
%                     with values
%                      'parameter', 'compartment', ''
%  timeVar:           Variable with time, assumed start of therapy is t=startTime.
%  startTime
%  myInterventionIDs: a cell array, 1 x n time points and variables
%
% RETURNS
% myBinTableRECIST

% TODO: Add proofing of inputs
continueFlag = true;

tableVariableNamesFixed = {'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID','PatientIDVar','TRTVar','BRSCOREVar','RSCOREVar'};
tableVariableNames = [tableVariableNamesFixed,{'weight','binEdge1','binEdge2','binEdge3','expN','expBin1','expBin2','expBin3','expBin4','predN','predBin1','predBin2','predBin3','predBin4'}];
myBinTableRECIST = cell2table(cell(0,length(tableVariableNames)));
myBinTableRECIST.Properties.VariableNames = tableVariableNames;

if continueFlag

    allExpDataIDs = getExpDataIDs(myWorksheet);

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
        % curData = curData(find(ismember(curData{:,BRSCOREVar},{'PD','SD','PR','CR'})),:);
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
                firstRow = find(curData{curRows,timeVar}==startTime);
                if isempty(firstRow)
                    firstRow = find(curData{curRows,timeVar} > startTime);
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
            % First get the bin edges from all variable data points
            allCurVarValues = selectData{curRows, curVar};
            myBinEdgeValues = myBinEdges{varCounter};
            if length(myBinEdgeValues) ~= 3
                myBinEdgeValues = [(median(allCurVarValues) - min(allCurVarValues))/2+min(allCurVarValues),median(allCurVarValues),(max(allCurVarValues) - median(allCurVarValues))/2+median(allCurVarValues)];
            end
            % Now we need to step through and evaluate for the times
            curTimes = unique(selectData{curRows,timeVar});
            for curTimeCounter = 1: length(curTimes)
                curDataRows = find(selectData{curRows,timeVar} == curTimes(curTimeCounter));
                if length(curDataRows) >= 2
                    curData = selectData{curRows(curDataRows),curVar};
                    % We only keep add a row if the largest value is bigger
                    % than the smallest.  Sometimes when normailizing to
                    % baseline variables with zero variability are created
                    % at T=0 and it makes sense to filter these.
                    if max(curData) > min(curData)
                        curProbs = wtdBinProb(curData', ones(1, length(curData))/length(curData), myBinEdgeValues);
                        expN = length(curData);
                        curRow = {curTimes(curTimeCounter), interventionID, mySimVars{varCounter}, mySimVarTypes{varCounter}, myExpDataID, timeVar, curVar, PatientIDVar, TRTVar, BRSCOREVar, RSCOREVar};
                        curRow = [curRow,{1,myBinEdgeValues(1), myBinEdgeValues(2), myBinEdgeValues(3), expN, curProbs(1), curProbs(2), curProbs(3), curProbs(4), nan, nan, nan, nan, nan}];
                        curRow = cell2table(curRow);
                        curRow.Properties.VariableNames = myBinTableRECIST.Properties.VariableNames; 
                        myBinTableRECIST = [myBinTableRECIST; curRow];
                    end
                end
            end
        end
    end
end
        
end