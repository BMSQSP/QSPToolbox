function myUpdatedObject = updateTablesWithSubpopformat(myInputObject)
% This is a function to update the tables attached to mapelOptions
% and VPop objects to convert to a new format of subpoptable.  
% The transition was made around in-house rev1228???.
%
% ARGUMENTS: 
%  myInputObject:     A VPop, VPopRECIST, VPopRECISTnoBin, 
%                      mapelOptions, mapelOptionsRECIST, or 
%                      mapelOptionsRECISTnoBin object to convert.
%
% RETURNS:
%  myUpdatedObject    A VPop, VPopRECIST, VPopRECISTnoBin, 
%                      mapelOptions, mapelOptionsRECIST, or 
%                      mapelOptionsRECISTnoBin object with
%                      the tables.

continueFlag = true;
if nargin > 1
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: a VPop, VPopRECIST, mapelOptions, or mapelOptionsRECIST.'])
    continueFlag = false;	
elseif nargin > 0
    continueFlag = true;
else
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: a VPop, VPopRECIST, mapelOptions, or mapelOptionsRECIST.'])
    continueFlag = false;
end

if ~(isa(myInputObject,'VPop') || isa(myInputObject,'VPopRECIST') || isa(myInputObject,'mapelOptions') || isa(myInputObject,'mapelOptionsRECIST') )
	warning(['Wrong input arguments for ',mfilename,'. Should provide: a VPop, VPopRECIST, mapelOptions, mapelOptionsRECIST.'])
	continueFlag = false;
end

myUpdatedObject = myInputObject;

if continueFlag
    updateFlag = false;
%     if isa(myInputObject, 'VPopRECIST') || isa(myInputObject, 'mapelOptionsRECIST')
%         updateTables = {'expData','mnSDTable','binTable','distTable','distTable2D','corTable','brTableRECIST','rTableRECIST'};
%     elseif isa(myInputObject, 'VPop') || isa(myInputObject, 'mapelOptions')
%         updateTables = {'expData','mnSDTable','binTable','distTable','distTable2D','corTable'};
%     end
    updateTables = myInputObject.subpopTable;
    [nRow, nCol] = size(updateTables);
    colnames = updateTables.Properties.VariableNames;
    for colCounter = 1:length(colnames)
        if ismember(colnames{colCounter},{'time','interventionID','elementID','elementType','comparator','value'})
            curCol = updateTables.(colnames{colCounter}); 
            updateCol = cell(nRow,1);
            if isequal(colnames{colCounter},'time') || isequal(colnames{colCounter},'value')
                if isnumeric(curCol)
                    updateCol{1} = num2cell(curCol(1));
                    updateFlag = true;
                end        
            end
            if nRow>1 
                if isequal(colnames{colCounter},'time') || isequal(colnames{colCounter},'value')
                    for rowi = 2:nRow
                        updateCol{rowi} = num2cell(curCol(rowi));
                    end
                else
                    for rowi = 2:nRow
                        updateCol(rowi) = num2cell(curCol(rowi));
                    end                    
                end
            end
            if updateFlag == true
                updateTables.(colnames{colCounter})=updateCol;  
            end
        end
    end
    myUpdatedObject.subpopTable =  updateTables;
    if updateFlag
        warning(['Detected old subpop format and updated ',class(myInputObject),' subpop tables to the new format.  Note that you will have to manually add a subpopTable to the object, which can be created by calling createSubpopTable with the appropriate worksheet.'])
    end
else
    warning(['Unable to run ',mfilename,'.  Returning input object.'])
end