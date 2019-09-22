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

optimizeType = myVPop.optimizeType;
% We don't adjust the parallel pool status for simplex
% since simplex cannot use parallel processing
% we may want to use the parallel pool to run multiple
% MAPEL runs in parallel.
if sum(ismember({'simplex'},optimizeType)) < 1
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end
	% We will use default pool settings
	mySimulateOptions = simulateOptions;
	mySimulateOptions = checkNWorkers(mySimulateOptions);
end

% Create the pool early, there are a few processes that
% may use it: both linearCalibraiton and the swarm optimization.
if sum(ismember({'simplex'},optimizeType)) < 1
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end    
    myPool = parpool(mySimulateOptions.clusterID,mySimulateOptions.nWorkers,'SpmdEnabled',false);
end

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
			initialPWs = getInitialPWs(myVPop, initialPWs, nAdd, mySimulateOptions.nWorkers*20);
		end
	end	
    
	[nTest, ~] = size(initialPWs);
	myPWTrans = nan(nTest,nTransPWs);
	for transCounter = 1 : nTest
		myPWTrans(transCounter,:)=hyperTransform(initialPWs(transCounter,:));
	end
end

if strcmp(myVPop.pwStrategy, 'bin')
    anonymousFunction = @(x)evaluateObjective(myVPop, x);
else
	anonymousFunction = @(x)evaluateObjectiveNoBin(myVPop, x);
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

% Clean up the worker pool
if sum(ismember({'simplex'},optimizeType)) < 1
	delete(myPool);
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