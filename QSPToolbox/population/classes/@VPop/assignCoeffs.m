function obj = assignCoeffs(obj, myWorksheet)
% Copy the VP coefficients from a worksheet to the `coeffsTable` property of the
% VPop object. 
% 
% EXAMPLE
% >> myVPop.assignCoeffs(myWorksheet)
%
% SEE ALSO
% VPop.coeffsTable
    mycoeffsTable = getVPCoeffs(myWorksheet);
    obj.coeffsTable = mycoeffsTable;
end
