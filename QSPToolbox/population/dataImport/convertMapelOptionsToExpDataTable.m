function myExpDataTable = convertMapelOptionsToExpDataTable(myWorksheet, myMapelOptions, patientIDVar)
% This function takes experimental data embedded in a worksheet and
% converts it to a data table format for use in MAPEL, guided by the
% data tables in the mapelOptions
%
% ARGUMENTS:
% myWorksheet:           a worksheet containing the data
% myMapelOptions
% patientIDVar:          variable with the patient IDs.  Currently only
%                        supported with RECIST calibrations
%
% RETURNS
% myExpDataTable

continueFlag = true;
if nargin > 3
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: myWorksheet, myMapelOptions, and patientIDVar if RECIST-based.'])
    continueFlag = false;
elseif ((nargin < 2) || ((nargin < 3)&&(isa(myMapelOptions,'mapelOptions'))))
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: myWorksheet, myMapelOptions, and patientIDVar if RECIST-based.'])
    continueFlag = false;
else
    continueFlag = true;    
end

if continueFlag
    if sum(ismember(class(myMapelOptions),{'mapelOptions','mapelOptionsRECIST','mapelOptionsRECISTnoBin'}))<0;
        warning(['Should provide: myWorksheet, myMapelOptions for ',mfilename,'.'])
        continueFlag = false;
    elseif sum(ismember(class(myMapelOptions),{'mapelOptions'}))<0
        warning(['The provided myMapelOptions for ',mfilename,' is not RECIST-based.  This function does not yet support 2D distributions with non-RECIST workflows.'])
    end
end

if continueFlag
    if isa(myMapelOptions, 'mapelOptions')
		%not yet supported
        %lookupRows  = cell2table(cell(0,7), 'VariableNames', {'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'PatientIDVar'});
		lookupRows  = cell2table(cell(0,6), 'VariableNames', {'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID'});
    else
        lookupRows  = cell2table(cell(0,10), 'VariableNames', {'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'PatientIDVar', 'TRTVar', 'BRSCOREVar', 'RSCOREVar'});
    end
    tablesToCheck = {'mnSDTable','binTable','distTable'};
    for tableCounter = 1: length(tablesToCheck)
        curTable = get(myMapelOptions,tablesToCheck{tableCounter});
        if isa(myMapelOptions, 'mapelOptions')
            %lookupRows = [lookupRows;curTable(:,{'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'PatientIDVar'})];
			lookupRows = [lookupRows;curTable(:,{'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID'})];
        else
            lookupRows = [lookupRows;curTable(:,{'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'PatientIDVar', 'TRTVar', 'BRSCOREVar', 'RSCOREVar'})];
        end
    end
    lookupRows = unique(lookupRows,'rows');
    [nRows, nCols] = size(lookupRows);
    for rowCounter = 1 : nRows
        if rowCounter == 1
            if isa(myMapelOptions, 'mapelOptions')
                %myExpDataTable = createExpDataTable(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1}, lookupRows{rowCounter,'PatientIDVar'}{1});
				myExpDataTable = createExpDataTable(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1});
            else
                myExpDataTable = createExpDataTableRECIST(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1}, lookupRows{rowCounter,'PatientIDVar'}{1}, lookupRows{rowCounter,'TRTVar'}, lookupRows{rowCounter,'BRSCOREVar'}{1}, lookupRows{rowCounter,'RSCOREVar'}{1});
            end
        else
            if isa(myMapelOptions, 'mapelOptions')
                %myExpDataTable = createExpDataTable(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1}, lookupRows{rowCounter,'PatientIDVar'}{1}, myExpDataTable);
				myExpDataTable = createExpDataTable(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1}, myExpDataTable);
            else
                myExpDataTable = createExpDataTableRECIST(myWorksheet, lookupRows{rowCounter,'interventionID'}{1}, lookupRows{rowCounter,'elementID'}{1}, lookupRows{rowCounter,'elementType'}{1}, lookupRows{rowCounter,'expDataID'}{1}, lookupRows{rowCounter,'expTimeVarID'}{1}, lookupRows{rowCounter,'expVarID'}{1}, lookupRows{rowCounter,'PatientIDVar'}{1}, lookupRows{rowCounter,'TRTVar'}, lookupRows{rowCounter,'BRSCOREVar'}{1}, lookupRows{rowCounter,'RSCOREVar'}{1}, myExpDataTable);
            end
        end    
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])
end    