function myVPop = initializeOptionPropertiesToVPop(myMapelOptions)
% This is a simple utility function to copy all of the properties from
% an options object to a VPop.
%
% ARGUMENTS
%  myMapelOptions:         (Required) an instance of a mapelOptions object.
%
% RETURNS
%  VPop:                    an instance of a VPop
%
if isa(myMapelOptions, 'mapelOptions')
    myVPop = VPop;
elseif isa(myMapelOptions, 'mapelOptionsRECIST')
    myVPop = VPopRECIST;
else
    myVPop = VPopRECISTnoBin;
end
% First, get the properties from myMapelOptions
% Copying over myMapelOptions.initialPWs is not done
% here.  It is assumed appropriate treatment of these
% is done in a calling script
myVPop.spreadOut = myMapelOptions.spreadOut;
myVPop.minIndPVal = myMapelOptions.minIndPVal;
myVPop.exactFlag = myMapelOptions.exactFlag;        
myVPop.useEffN = myMapelOptions.useEffN;    
myVPop.optimizeTimeLimit = myMapelOptions.optimizeTimeLimit;
myVPop.optimizeType = myMapelOptions.optimizeType;        
myVPop.optimizePopSize = myMapelOptions.optimizePopSize;
myVPop.objectiveLimit = myMapelOptions.objectiveLimit; 	
myVPop.intSeed = myMapelOptions.intSeed;
myVPop.nIters = myMapelOptions.nIters;
myVPop.tol = myMapelOptions.tol;
myVPop.expData = myMapelOptions.expData;
myVPop.mnSDTable = myMapelOptions.mnSDTable;
myVPop.binTable = myMapelOptions.binTable;
myVPop.distTable = myMapelOptions.distTable; 
myVPop.distTable2D = myMapelOptions.distTable2D;
myVPop.corTable = myMapelOptions.corTable;
myVPop.minEffN = myMapelOptions.minEffN;

if isa(myVPop,'VPopRECIST')

		myVPop.brTableRECIST = myMapelOptions.brTableRECIST;
		myVPop.rTableRECIST = myMapelOptions.rTableRECIST;        
        myVPop.relSLDvar = myMapelOptions.relSLDvar;
        myVPop.absALDVar = myMapelOptions.absALDVar;
        myVPop.crCutoff = myMapelOptions.crCutoff;        
end	    
if isa(myVPop,'VPopRECISTnoBin')
        myVPop.distTable2D = myMapelOptions.distTable2D;
		myVPop.corTable = myMapelOptions.corTable;
		myVPop.brTableRECIST = myMapelOptions.brTableRECIST;
		myVPop.rTableRECIST = myMapelOptions.rTableRECIST;        
        myVPop.relSLDvar = myMapelOptions.relSLDvar;
        myVPop.absALDVar = myMapelOptions.absALDVar;
        myVPop.crCutoff = myMapelOptions.crCutoff;        
end
end