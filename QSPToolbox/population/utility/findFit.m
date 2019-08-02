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

if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPop')
	initialBinProbs = myVPop.binProbs;
	[myNAxis, myNBins] = size(initialBinProbs);
	myTransProbs = nan(myNAxis, myNBins-1);
	for axisCounter = 1 : myNAxis
		myTransProbs(axisCounter,:) = hyperTransform(initialBinProbs(axisCounter,:));
	end
	transProbVect = probsToProbVect(myTransProbs);
	nTransProbs = length(transProbVect);
elseif isa(myVPop,'VPopRECISTnoBin')
	initialPWs=myVPop.pws;
	nTransPWs=length(initialPWs)-1;
	[nTest, ~] = size(initialPWs);
	
	if ~(ismember({'simplex'},optimizeType))
		% We will also try supplementing the VP scores scaled into PWs.  This
		% likely won't be very effective at finding optimal solutions but will add
		% points where VPs in sparser regions of the distributions relative to the data 
		% are more highly weighted which could help to weight to where it is needed
		% without overly focusing on a few VPs
		vpScores = scoreWorksheetVPs(myVPop,1:(nTransPWs+1),1:(nTransPWs+1));
		vpScores = (vpScores + 1E-12)./sum((vpScores + 1E-12),2);
		
        % Get an initial point from linear calibrate.  First use a uniform
        % assumption with bagging, optimize to get a prior
        % prevalence weight assumption then re-optimize with bagging.
        myOptimOptions = LinearCalibrationOptions;
        myOptimOptions.cdfProbsToFit = 0.05:0.05:0.95;
        myOptimOptions.optimizationAlgorithm = "nnls";
		myOptimOptions.optimizationAlgorithmOptions.Accy = 0;
        myOptimOptions.priorPrevalenceWeightAssumption = "uniform";
        myOptimOptions.nBootstrapIterations = mySimulateOptions.nWorkers*20;
        myOptimOptions.fractionVPsPerBaggingIteration=.5;
        myOptimOptions.method = "bagging";
        linearCalibrationObject = LinearCalibration(myVPop,'optimOptions',myOptimOptions);
        try
            linearCalibrationObject = linearCalibrationObject.run();
            if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
                linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
            else
                linearCalibrationPWs = '';
            end
        catch
            % This will sometimes fail if all of the bagged
            % trials don't run.
            linearCalibrationPWs = '';
        end
        if isnumeric(linearCalibrationPWs)
            myOptimOptions.priorPrevalenceWeightAssumption = "specified";		
            myOptimOptions.nBootstrapIterations = mySimulateOptions.nWorkers*20;	
            myOptimOptions.fractionVPsPerBaggingIteration=.5;
            myOptimOptions.method = "bagging";
            linearCalibrationObject = LinearCalibration(linearCalibrationObject.OptimizedVPop,'optimOptions',myOptimOptions);
            try
                linearCalibrationObject = linearCalibrationObject.run();
                if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
                    linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
                else
                    % Dummy statement, keep the original
                    linearCalibrationPWs = linearCalibrationPWs;
                end
            catch
                % This will sometimes fail if all of the bagged 
                % trials don't run.
                % Dummy statement, keep the original
                linearCalibrationPWs = linearCalibrationPWs;
            end
        end
        
        % We reserve about half of the allowed particles
		% for fully random assignment
		initialPWs = [initialPWs;vpScores];
		if isnumeric(linearCalibrationPWs)
			initialPWs = [initialPWs;linearCalibrationPWs];
            [nTest, ~] = size(initialPWs);
            nAddOld = max(floor(myVPop.optimizePopSize/4-nTest/2),0);
            nAddLin = nAddOld;
        else
            [nTest, ~] = size(initialPWs);
            nAddOld = floor(myVPop.optimizePopSize/2-nTest);
            nAddLin = 0;
        end
		
		[nTest, ~] = size(initialPWs);
        % Proof again to make sure we don't exceed
		% the allowed swarm size.
		nTest = min(nTest,myVPop.optimizePopSize);
		initialPWs = initialPWs(1:nTest,:);
		% Try supplementing with initial guesses that weight VPs near
		% the initial points
		if (((myVPop.optimizePopSize-nTest-nAddOld) > 0) && (nAddOld>0))
			vpScores = spreadPWsToNeighbors(initialPWs(1,:), myVPop.coeffsTable, nAddOld);
			initialPWs = [initialPWs;vpScores(2:end,:)];
        end
        % Also try supplementing with PW solutions near the bagged linear
        % calibrate optimal solution
        if (((myVPop.optimizePopSize-nTest-nAddLin) > 0) && (nAddLin>0))
			vpScores = spreadPWsToNeighbors(linearCalibrationObject.OptimizedVPop.pws, myVPop.coeffsTable, nAddLin);
			initialPWs = [initialPWs;vpScores(2:end,:)];
        end      
        % Update again based on the additions.
        [nTest, ~] = size(initialPWs);
	end	
    
	myPWTrans = nan(nTest,nTransPWs);
	for transCounter = 1 : nTest
		myPWTrans(transCounter,:)=hyperTransform(initialPWs(transCounter,:));
	end
end

if isa(myVPop,'VPop')
    anonymousFunction = @(x)evaluateObjective(myVPop, x);
elseif isa(myVPop,'VPopRECIST')
    anonymousFunction = @(x)evaluateObjectiveRECIST(myVPop, x);
else
	anonymousFunction = @(x)evaluateObjectiveRECISTnoBin(myVPop, x);
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
	if ~isa(myVPop,'VPopRECISTnoBin')
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
	if ~isa(myVPop,'VPopRECISTnoBin')
		optimOptions.PopInitRange = cat(1,ones(1,nTransProbs)*-pi/2,ones(1,nTransProbs)*pi/2);
		optimOptions.InitialPopulation = transProbVect;
		optimOptions.MutationFcn = {@mutationuniform, 1.5/nTransProbs};
		[optTransProbVect,fVal,exitFlag,output] = ga(anonymousFunction,nTransProbs,[],[],[],[],ones(1,nTransProbs)*-pi/2,ones(1,nTransProbs)*pi/2,[],optimOptions);
		
	else
		optimOptions.PopInitRange = cat(1,ones(1,nTransPWs)*-pi/2,ones(1,nTransPWs)*pi/2);
		optimOptions.InitialPopulation = myPWTrans;
		optimOptions.MutationFcn = {@mutationuniform, 1.5/nTransPWs};	
		[optTransPWsVect,fVal,exitFlag,output] = ga(anonymousFunction,nTransPWs,[],[],[],[],ones(1,nTransPWs)*-pi/2,ones(1,nTransPWs)*pi/2,[],optimOptions);

	end	
end

% Clean up the worker pool
if sum(ismember({'simplex'},optimizeType)) < 1
	delete(myPool);
end
if ~isa(myVPop,'VPopRECISTnoBin')     
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