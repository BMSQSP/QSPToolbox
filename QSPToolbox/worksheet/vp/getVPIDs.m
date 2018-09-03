function vpIDs = getVPIDs(myWorksheet)
% Get VP IDs from a model
% ARGUMENTS
% myWorksheet: a worksheet
%
% RETURNS
% vpIDs: an array of cells with the vp names

VPs = myWorksheet.vpDef;
vpIDs = cell(1, length(VPs));
for theVPcounter = 1 : length(VPs)
    curVP = VPs{theVPcounter};
    vpIDs{theVPcounter} = curVP.('ID');
end
        
        
end