function myVPop = findFit(myVPop)
% This is main function for optimizing the fit to data.
% This function operates on and returns a VPop object
% but is a separate function rather than a VPop method.
%
% ARGUMENTS
%  myVPop: An object of class VPop, VPopRECIST, or VPopRECISTnoBin
%
% RETURNS
%  myVPop: this function will populate several properties:
%          binProbs (if not a VPopRECISTnoBin)
%          pws
%          mnSDTable (add VPop predictions)
%          binTable (add VPop predictions)
%          distTable (add VPop predictions)
%          gofMn
%          gofSD
%          gofBin
%          gofDist
%          gof
%

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolSize = 1;
else
    poolSize = p.NumWorkers;
end

optimizeType = myVPop.optimizeType;
if strcmp(myVPop.pwStrategy, 'bin')
	initialBinProbs = myVPop.binProbs;
	[myNAxis, myNBins] = size(initialBinProbs);
	myTransProbs = nan(myNAxis, myNBins-1);
	for axisCounter = 1 : myNAxis
		myTransProbs(axisCounter,:) = hyperTransform(initialBinProbs(axisCounter,:));
	end
	transProbVect = probsToProbVect(myTransProbs);
	nTransProbs = length(transProbVect);
else
	initialPWs=myVPop.pws;
	nTransPWs=length(initialPWs)-1;
	
	if ~(ismember({'simplex'},optimizeType))
		[nTest, ~] = size(initialPWs);
		nAdd = floor(myVPop.optimizePopSize/2-nTest);
		if nAdd > 0
            nBoostrapIterations = max(100,poolSize*5);
			initialPWs = getInitialPWs(myVPop, initialPWs, nAdd, nBoostrapIterations);
		end
	end	
    
	[nTest, ~] = size(initialPWs);
	myPWTrans = nan(nTest,nTransPWs);
	for transCounter = 1 : nTest
		myPWTrans(transCounter,:)=hyperTransform(initialPWs(transCounter,:));
	end
end

% We'll strip away unnecessary data for the workers
% before running the swarm optimization.
% They don't need:
% expData
% RECISTSimFilters
% since the tables already have extracted the needed exp data
% The tables are updated with simData during the iterations
% but the sim filters can at least be dropped at this point
% if present.
% Stripping these fields were not not noted
% to make much of a difference in performance,
% but this is done as a precautionary measure.
computeVPop=myVPop;
computeVPop.expData=[];
if isa(myVPop,'VPopRECIST')
    computeVPop.recistSimFilter=[];
end

if strcmp(myVPop.pwStrategy, 'bin')
    anonymousFunction = @(x)evaluateObjective(computeVPop, x);
else
	anonymousFunction = @(x)evaluateObjectiveNoBin(computeVPop, x);
end

if sum(ismember({'simplex'},optimizeType)) > 0
    optimOptions = optimset('fminsearch');    
    optimOptions.Display = 'iter';    
    optimOptions.MaxIter = myVPop.nIters;
    optimOptions.MaxFunEvals = myVPop.nIters;    
    % Parallel is not supported for fminsearch
    optimOptions.UseParallel = false;
    optimOptions.TolFun = myVPop.tol;
    optimOptions.TolX = myVPop.tol;
	optimOptions.ObjectiveLimit = myVPop.objectiveLimit;
	if ~isa(myVPop,'VPopRECISTnoBin')
		[optTransProbVect,fVal,exitFlag,output] = fminsearch(anonymousFunction,transProbVect,optimOptions);
	else
		[optTransPWsVect,fVal,exitFlag,output] = fminsearch(anonymousFunction,myPWTrans,optimOptions);
	end
elseif sum(ismember({'pso'},optimizeType)) > 0
    optimOptions = optimoptions('particleswarm');
    optimOptions.Display = 'iter';
    optimOptions.MaxTime = myVPop.optimizeTimeLimit;
    optimOptions.TolFun = myVPop.tol;
    optimOptions.UseParallel = true;                
    optimOptions.SwarmSize = myVPop.optimizePopSize;  
	optimOptions.ObjectiveLimit = myVPop.objectiveLimit;
	if strcmp(myVPop.pwStrategy, 'bin')	
		optimOptions.InitialSwarm = transProbVect;
		[optTransProbVect,fVal,exitFlag,output] = particleswarm(anonymousFunction,nTransProbs,ones(1,nTransProbs)*-pi/2,ones(1,nTransProbs)*pi/2,optimOptions);
	else
		optimOptions.InitialSwarm = myPWTrans;
        [optTransPWsVect,fVal,exitFlag,output] = particleswarm(anonymousFunction,nTransPWs,ones(1,nTransPWs)*-pi/2,ones(1,nTransPWs)*pi/2,optimOptions); 		
		
	end
else
    optimOptions = gaoptimset;
    optimOptions.Display = 'diagnose';
    optimOptions.MaxTime = myVPop.optimizeTimeLimit;
    optimOptions.TolFun = myVPop.tol;
    optimOptions.UseParallel = true;    
    optimOptions.PopulationSize = myVPop.optimizePopSize;
	optimOptions.ObjectiveLimit = myVPop.objectiveLimit;
	if strcmp(myVPop.pwStrategy, 'bin')	
		optimOptions.PopInitRange = cat(1,ones(1,nTransPWs)*-pi/2,ones(1,nTransPWs)*pi/2);
		optimOptions.InitialPopulation = myPWTrans;
		optimOptions.MutationFcn = {@mutationuniform, 1.5/nTransPWs};	
		[optTransPWsVect,fVal,exitFlag,output] = ga(anonymousFunction,nTransPWs,[],[],[],[],ones(1,nTransPWs)*-pi/2,ones(1,nTransPWs)*pi/2,[],optimOptions);		
	else
		optimOptions.PopInitRange = cat(1,ones(1,nTransProbs)*-pi/2,ones(1,nTransProbs)*pi/2);
		optimOptions.InitialPopulation = transProbVect;
		optimOptions.MutationFcn = {@mutationuniform, 1.5/nTransProbs};
		[optTransProbVect,fVal,exitFlag,output] = ga(anonymousFunction,nTransProbs,[],[],[],[],ones(1,nTransProbs)*-pi/2,ones(1,nTransProbs)*pi/2,[],optimOptions);
	end	
end

if strcmp(myVPop.pwStrategy, 'bin')	    
	myProbTrans = transpose(reshape(optTransProbVect, myNBins-1, myNAxis));
	myBinProbs = nan(myNAxis, myNBins);
	for axisCounter = 1 : myNAxis
		myBinProbs(axisCounter,:) = invHyperTransform(myProbTrans(axisCounter,:));
	end
	myVPop.binProbs = myBinProbs;
	myVPop = myVPop.assignPWs();
else
	myVPop.pws=invHyperTransform(optTransPWsVect);
end

% Next update table predicted values
myVPop = myVPop.addPredTableVals();
% Next update the individual GOF statistics
myVPop = evaluateGOF(myVPop);
end    