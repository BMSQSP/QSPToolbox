function myVPop = restartMapel(myVPop, myRandomStart)
% This is a function to start developing a VPop using data in an existing
% existing VPop as a starting point.
%
% ARGUMENTS
%  myVPop:           A VPop instance.  Note restarting will fail
%                    if there is not enough information in the VPop
%                    from a previous run.
%  myRandomStart:    (optional) A non-negative numeric value. If 0, the initial 
%                    parameters are not perturbed in anyway.  Otherwise the 
%                    transformed initial bin probabilities are perturbed  
%                    with random (normally distributed) noise that has
%                    normalized standard deviation of the specified value.  
%                    If no initialProbs are specified, a value of zero 
%                    results in uniform bin probabilities and a value >0 
%                    results in uniform normal initial bin probabilities  
%                    (renormalized so the margin is 1).  Default is 0.
%
% Returns
%  VPop:             an instance of a VPop
%
continueFlag = false;
if nargin > 2
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: myVpop and optionally myRandomStart.'])
    continueFlag = false;
elseif nargin > 1
    continueFlag = true;
elseif nargin > 0
    myRandomStart = 0;
    continueFlag = true;
else
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: myVpop and optionally myRandomStart.'])
end

if continueFlag
    if (~strcmp(class(myVPop),'VPop') & ~strcmp(class(myVPop),'VPopRECIST') & ~strcmp(class(myVPop),'VPopRECISTnoBin'))
        continueFlag = false;
        warning(['Input VPop not recongized in call to ',mfilename,'.  Requires: myVpop and optionally myRandomStart.'])
    end
    if ~(myRandomStart>=0)
        continueFlag = false;
        warning(['Input myRandomStart not recongized in call to ',mfilename,'.  myRandomStart must be a non-negative number.'])
    end        
end

if continueFlag
    % We'll just call MAPEL, so all we need to do is feed the input VPop
    % back to a mapelOptions object.
	if strcmp(class(myVPop),'VPopRECIST')
		myMapelOptions = mapelOptionsRECIST;
	elseif strcmp(class(myVPop),'VPop')
		myMapelOptions = mapelOptions;
	elseif strcmp(class(myVPop),'VPopRECISTnoBin')
		myMapelOptions = mapelOptionsRECISTnoBin;
	end
    myMapelOptions.expData = myVPop.expData;
    myMapelOptions.mnSDTable = myVPop.mnSDTable;
    myMapelOptions.binTable = myVPop.binTable;
    myMapelOptions.distTable = myVPop.distTable;	
    myMapelOptions.optimizeType = myVPop.optimizeType;
    myMapelOptions.optimizeTimeLimit = myVPop.optimizeTimeLimit;
    myMapelOptions.nIters = myVPop.nIters;
    if ~strcmp(class(myVPop),'VPopRECISTnoBin')
    [nAxis, nBins] = size(myVPop.binProbs);
        myMapelOptions.nBins = nBins;
        myMapelOptions.initialProbs = myVPop.binProbs;
    end
    if strcmp(class(myVPop),'VPopRECISTnoBin')
        myMapelOptions.initialPWs = myVPop.pws;
    end    
    myMapelOptions.tol = myVPop.tol;
    myMapelOptions.spreadOut = myVPop.spreadOut;
    myMapelOptions.minIndPVal = myVPop.minIndPVal;	
    myMapelOptions.exactFlag = myVPop.exactFlag;    
    myMapelOptions.minEffN = myVPop.minEffN;
    myMapelOptions.optimizePopSize = myVPop.optimizePopSize;  
	myMapelOptions.objectiveLimit = myVPop.objectiveLimit; 
    myMapelOptions.optimizeTimeLimit = myVPop.optimizeTimeLimit;    
    myMapelOptions.randomStart = myRandomStart;
    myMapelOptions.intSeed = myVPop.intSeed;
	myMapelOptions.useEffN = myVPop.useEffN;
    myMapelOptions.distTable2D = myVPop.distTable2D;
	myMapelOptions.corTable = myVPop.corTable;
	if strcmp(class(myVPop),'VPopRECIST') | strcmp(class(myVPop),'VPopRECISTnoBin')		
		myMapelOptions.brTableRECIST = myVPop.brTableRECIST;
		myMapelOptions.rTableRECIST = myVPop.rTableRECIST;        
        myMapelOptions.relSLDvar = myVPop.relSLDvar;
        myMapelOptions.absALDVar = myVPop.absALDVar;
        myMapelOptions.crCutoff = myVPop.crCutoff;        
    end
	myVPop = mapel(myVPop, myMapelOptions);
else
    warning(['Unable to run ',mfilename,'.  Returning myVPop.'])
end
end