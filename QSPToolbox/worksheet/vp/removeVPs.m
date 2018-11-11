function myWorksheet = removeVPs(myWorksheet, vpIDs)
% Remove VPs from a worksheet based on ID
% ARGUMENTS
% myWorksheet: a worksheet
% vpIDs: ids of VPs to remove from worksheet
%
% RETURNS
% myWorksheet: an updated worksheet without the VPs

allVPids = getVPIDs(myWorksheet);
nOriginalVPs = length(allVPids);
curIndices = (find(ismember(allVPids,vpIDs)));
if length(curIndices) > 0
    [size1, size2] = size(myWorksheet.results);
    if size2 == nOriginalVPs;
        myWorksheet.results(:,curIndices) = [];
    end
    [size1, size2] = size(myWorksheet.axisProps.axisVP.coefficients);
    if size2 == nOriginalVPs;
        myWorksheet.axisProps.axisVP.coefficients(:,curIndices) = [];
    end
    myWorksheet.vpDef(curIndices) = [];    
end
    
end
