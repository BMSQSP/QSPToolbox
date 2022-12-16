function myDistTable = createDistTable(myWorksheet,myExpDataIDs,PatientIDVar,myVars,mySimVars,mySimVarTypes,timeVar,startTime,myInterventionIDs)
% This function takes experimental data embedded in a worksheet and
% converts it to a distribution table for nonRECIST mapelOptions, bypassing reading from a mapelOptions
% "experimental data table" and directly taking data from the worksheet
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
%  myDistTable

continueFlag = true;

commonNames = loadCommonNames();
tableVariableNamesFixed = commonNames.VPOPTABLEVARNAMESFIXED;
tableVariableNames = [tableVariableNamesFixed,{'weight','expN', 'expSample', 'predN', 'predIndices','predSample', 'predProbs','expCombinedIndices','simCombinedIndices','combinedPoints'}];
myDistTable = cell2table(cell(0,length(tableVariableNames)));
myDistTable.Properties.VariableNames = tableVariableNames;

% Proofing of the input
% variables starts here
if continueFlag
	testN = length(myExpDataIDs);
	if sum([length(myExpDataIDs),length(myVars),length(mySimVars),length(mySimVarTypes),length(myInterventionIDs)] ~= testN) > 0
		continueFlag = false;
		warning(['Not all input cell array arguments are of consistent length in ',mfilename,'.  Exiting.'])
	end
end

if continueFlag
    
    allExpDataIDs = getExpDataIDs(myWorksheet);
    allInterventionIDs = getInterventionIDs(myWorksheet);
    
    for varCounter = 1 : testN
    
        curVar = myVars{varCounter};
        myExpDataID = myExpDataIDs{varCounter};
        interventionID = myInterventionIDs{varCounter};
  
        curIndex = find(ismember(allExpDataIDs,myExpDataID));
        
        if isempty(curIndex)
            warning(['Not able to find experimental dataset ',myExpDataID,' in ',mfilename,'.'])
            continueFlag = false;
        else        
            curData = myWorksheet.expData{curIndex}.data;  % nonRECIST: no off-treatment info, take all data
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
	
	% Also proof whether the intervention/model variable inputs are unique.
	% Non-uniqueness would suggest mapping the same model/intervention
	% combination onto multiple outputs, which should not be allowed	
	mySimVars = reshape(mySimVars,testN,1);
	mySimVarTypes = reshape(mySimVarTypes,testN,1);
	myInterventionIDs = reshape(myInterventionIDs,testN,1);
    myExpDataIDs = reshape(myExpDataIDs,testN,1);
	combineRows = [mySimVars,mySimVarTypes,myInterventionIDs,myExpDataIDs];
	combineRows = cell2table(combineRows);
	[nRows,~] = size(unique(combineRows,'rows'));
	if nRows ~= testN
		warning(['Not all model variables, variable types, intervention sets are unique in call to ',mfilename,'.  Exiting.'])
		continueFlag = false;    
	end
end
     

if continueFlag
    for varCounter = 1 : testN
        curVar = myVars{varCounter};
        myExpDataID = myExpDataIDs{varCounter};
        interventionID = myInterventionIDs{varCounter};
        curIndex = find(ismember(allExpDataIDs,myExpDataID));
        curData = myWorksheet.expData{curIndex}.data;        
            
        selectData = curData;    % no filter based on RECIST, take all data
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
						% When initially transferring the data over,
						% assume it is from the "all" subpop.						
                        newRow = {1, curTimes(curTimeCounter), interventionID, mySimVars{varCounter}, mySimVarTypes{varCounter}, myExpDataID, timeVar, curVar, PatientIDVar};
                        newRow = [newRow,{1, length(curData), {sort(curData,'ascend')}, nan, {nan}, {nan}, {nan},{nan},{nan},{nan}}];
                        newRow = cell2table(newRow);
                        %newRow = table(curTimes(curTimeCounter), interventionID, mySimVars{varCounter}, mySimVarTypes{varCounter}, myExpDataID, timeVar, curVar,PatientIDVar,TRTVar,BRSCOREVar,RSCOREVar1, length(curData), {sort(curData,'ascend')}, nan, {nan}, {nan}, {nan},{nan},{nan},{nan})
                        newRow.Properties.VariableNames = myDistTable.Properties.VariableNames; 
                        myDistTable = [myDistTable; newRow];
                    end
                end
            end
        end
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])    
end    