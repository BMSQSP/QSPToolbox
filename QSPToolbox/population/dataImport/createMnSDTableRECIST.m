function myMnSDTableRECIST = createMnSDTableRECIST(myWorksheet,myExpDataIDs,PatientIDVar,TRTVar,BRSCOREVar,RSCOREVar,myVars,mySimVars,mySimVarTypes,timeVar,startTime,myInterventionIDs)
% This function takes experimental data embedded in a worksheet and
% converts it to a mean/sd table, bypassing the "experimental data table"
% directly accounting for whether patients should be excluded because of
% CR, PD, or they are clearly of treatment.
%
% ARGUMENTS:
%  myWorksheet:       a workshet with the experimental data
%  myExpDataIDs:      a cell array, 1 x n time points and variables
%                      with the experimental data IDs
%  PatientIDVar:      variable with the patient IDs
%  TRTVar:            a binary variable indicating whether a patient is on
%                      treatment
%  BRSCOREVar:        best overall response to date variable.
%  RSCOREVar:         current response
%  myVars:            a cell array, 1 x n time points and variables
%                      variables to get mean/SD data for.
%  mySimVars:         simulation variable names, map to myVars
%  mySimVarTypes:     a cell array, 1 x n time points and variables
%                     with values
%                      'parameter', 'compartment', ''
%  timeVar:           Variable with time, 
%                     NOTE: it is assumed start of therapy is t=startTime.
%  startTime
%  myInterventionIDs: a cell array, 1 x n time points and variables
%
% RETURNS
% myMnSDTableRECIST

% TODO: Add proofing of inputs
continueFlag = true;

tableVariableNamesFixed = {'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID','PatientIDVar','TRTVar','BRSCOREVar','RSCOREVar'};
tableVariableNames = [tableVariableNamesFixed,{'weightMean', 'weightSD', 'expN', 'expMean', 'expSD', 'predN', 'predMean', 'predSD'}];
myMnSDTableRECIST = cell2table(cell(0,length(tableVariableNames)));
myMnSDTableRECIST.Properties.VariableNames = tableVariableNames;

if continueFlag
    
    allExpDataIDs = getExpDataIDs(myWorksheet);
    expN = nan(0,1);
    expMean = nan(0,1);
    expSD = nan(0,1);
    expTime = nan(0,1);
    interverventionIDCell = cell(0,1);
    elementIDCell = cell(0,1);
    elementTypeCell = cell(0,1);
    expDataIDCell = cell(0,1);
    expTimeVarIDCell = cell(0,1);
    expVarIDCell = cell(0,1);  
    PatientIDCell = cell(0,1);
    TRTCell = cell(0,1);
    BRSCORECell = cell(0,1);
    RSCORECell = cell(0,1);
    
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
        % We want to scan down patient columns and keep data from patients
        % that have valid RECIST measures.
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

        % This code was obseleted by the previous lines
        % We will select response rows
        %responseRows = find(ismember(curData{:,BRSCOREVar},{'PD','SD','PR','CR'}));
        % We will also keep baseline rows, where response may not be measured
        % for now assume timeVar == 0 here
        %timeRows = find(curData{:,timeVar}==0);
        %responseRows = sort(unique([responseRows;timeRows]),'ascend');
        % curData = curData(responseRows,:);
        % Now also filter once a patient hits CR, PD

        % Now we create mean and SD table from data in selectData
        % Table will have length sum(nTimePointsVari)


        curRows = find(~isnan(selectData{:,curVar}));
        if length(curRows) >= 2
            curTimes = unique(selectData{curRows,timeVar});
            for curTimeCounter = 1: length(curTimes)
                curDataRows = find(selectData{curRows,timeVar} == curTimes(curTimeCounter));
                curData = selectData{curRows(curDataRows),curVar};
                if length(curDataRows) >= 2
                    expN = [expN; length(curDataRows)];
                    expMean = [expMean; mean(curData)];
                    expSD = [expSD; std(curData)];
                    expTime = [expTime; curTimes(curTimeCounter)];
                    interverventionIDCell = [interverventionIDCell; interventionID];
                    elementIDCell = [elementIDCell; mySimVars{varCounter}];
                    elementTypeCell = [elementTypeCell; mySimVarTypes{varCounter}];
                    expDataIDCell = [expDataIDCell; myExpDataID];
                    expTimeVarIDCell = [expTimeVarIDCell; timeVar];
                    expVarIDCell = [expVarIDCell; curVar];
                    PatientIDCell = [PatientIDCell; PatientIDVar];
                    TRTCell = [TRTCell; TRTVar];
                    BRSCORECell = [BRSCORECell; BRSCOREVar];
                    RSCORECell = [RSCORECell; RSCOREVar];                    
                end
            end
        end
    end
    
    myMnSDTableRECIST = table(expTime,interverventionIDCell,elementIDCell,elementTypeCell,expDataIDCell,expTimeVarIDCell,expVarIDCell,PatientIDCell,TRTCell,BRSCORECell,RSCORECell,ones(length(expTime),1),ones(length(expTime),1),expN,expMean,expSD,nan(length(expTime),1),nan(length(expTime),1),nan(length(expTime),1));
    % We also filter where SD > 0.  SD = 0 can show up when normalizing
    % to a baseline value, and it doesn't make sense to calibrate
    % the baseline value to zero.
    myMnSDTableRECIST.Properties.VariableNames = tableVariableNames;
    myMnSDTableRECIST = myMnSDTableRECIST(find(myMnSDTableRECIST{:,'expSD'}>0),:);
    end
        
end    