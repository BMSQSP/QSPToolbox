function initialPWs = getInitialPWs(myVPop,initialPWs, nAdd, nBootstrapIterations)
% This is a utility function to find additional starting
% points for the hypothesis-test driven swarming approach.
%
% ARGUMENTS
%  myVPop:               An object of VPopRECISTnoBin.  Other 
%                         VPop classes not yet supported.
%  initialPWs:           Starting initial PWs.  
%  nAdd:                 Number of PWs to add.
%  nBootstrapIterations: For linear calibration
% RETURNS
%  initialPWs:           An updated set of initial prevalence weights.
%

% We will also try supplementing the VP scores scaled into PWs.  This
% likely won't be very effective at finding optimal solutions but will add
% points where VPs in sparser regions of the distributions relative to the data 
% are more highly weighted which could help to weight to where it is needed
% without overly focusing on a few VPs
[~,nTransPWs] = size(myVPop.pws);
vpScores = scoreWorksheetVPs(myVPop,1:nTransPWs,1:nTransPWs);
vpScores = (vpScores + 1E-12)./sum((vpScores + 1E-12),2);

[nAddScores, ~] = size(vpScores);
nAddScores = min(nAddScores,nAdd);
nAdd = nAdd - nAddScores;

% Hold on to VP scores until later.
% Iteratively run linear calibrations
% to get more starting points
% with alternate assumptions for
% simulated data weight
testVPop=myVPop;
if nAdd >= 1
   myOptimOptions = LinearCalibrationOptions();
   myOptimOptions.cdfProbsToFit = 0.05:0.05:0.95;
   myOptimOptions.optimizationAlgorithm = "nnls";
   myOptimOptions.optimizationAlgorithmOptions.Accy = 0;
   myOptimOptions.priorPrevalenceWeightAssumption = "uniform";
   myOptimOptions.nBootstrapIterations = nBootstrapIterations;
   myOptimOptions.fractionVPsPerBaggingIteration=.5;
   myOptimOptions.method = "bagging";
   simN = 1;
   myOptimOptions.expWeightFuncHandle = @(expN, expSTD, dataGroupDescription) calculateExpWeight(expN, expSTD, dataGroupDescription, simN);
   linearCalibrationObject = LinearCalibration(testVPop,'optimOptions',myOptimOptions);
   try
      linearCalibrationObject = linearCalibrationObject.run('closeParallelPoolWhenFinished',false);
      if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
		 testVPop = linearCalibrationObject.OptimizedVPop;
		 myOptimOptions.priorPrevalenceWeightAssumption = "specified";
         linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
		 initialPWs = [initialPWs;linearCalibrationPWs];
		 nAdd = nAdd - 1;
       else
         linearCalibrationPWs = '';
       end
   catch
       % This will sometimes fail if all of the bagged
       % trials don't run.
       linearCalibrationPWs = '';
   end	   
end
if nAdd >= 1
   simN = 10;
   myOptimOptions.expWeightFuncHandle = @(expN, expSTD, dataGroupDescription) calculateExpWeight(expN, expSTD, dataGroupDescription, simN);
   linearCalibrationObject = LinearCalibration(testVPop,'optimOptions',myOptimOptions);
   try
      linearCalibrationObject = linearCalibrationObject.run('closeParallelPoolWhenFinished',false);
      if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
		 testVPop = linearCalibrationObject.OptimizedVPop;
		 myOptimOptions.priorPrevalenceWeightAssumption = "specified";
         linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
		 initialPWs = [initialPWs;linearCalibrationPWs];
		 nAdd = nAdd - 1;
       else
         linearCalibrationPWs = '';
       end
   catch
       % This will sometimes fail if all of the bagged
       % trials don't run.
       linearCalibrationPWs = '';
   end	   
end
if nAdd >= 1
   simN = 100;
   myOptimOptions.expWeightFuncHandle = @(expN, expSTD, dataGroupDescription) calculateExpWeight(expN, expSTD, dataGroupDescription, simN);
   linearCalibrationObject = LinearCalibration(testVPop,'optimOptions',myOptimOptions);
   try
      linearCalibrationObject = linearCalibrationObject.run('closeParallelPoolWhenFinished',false);
      if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
		 testVPop = linearCalibrationObject.OptimizedVPop;
		 myOptimOptions.priorPrevalenceWeightAssumption = "specified";
         linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
		 initialPWs = [initialPWs;linearCalibrationPWs];
		 nAdd = nAdd - 1;
       else
         linearCalibrationPWs = '';
       end
   catch
       % This will sometimes fail if all of the bagged
       % trials don't run.
       linearCalibrationPWs = '';
   end	   
end
if nAdd >= 1
   simN = 1000;
   myOptimOptions.expWeightFuncHandle = @(expN, expSTD, dataGroupDescription) calculateExpWeight(expN, expSTD, dataGroupDescription, simN);
   linearCalibrationObject = LinearCalibration(testVPop,'optimOptions',myOptimOptions);
   try
      linearCalibrationObject = linearCalibrationObject.run('closeParallelPoolWhenFinished',false);
      if isnumeric(linearCalibrationObject.OptimizedVPop.pws) && (min(linearCalibrationObject.OptimizedVPop.pws)>=0)
		 testVPop = linearCalibrationObject.OptimizedVPop;
		 myOptimOptions.priorPrevalenceWeightAssumption = "specified";
         linearCalibrationPWs = linearCalibrationObject.OptimizedVPop.pws;
		 initialPWs = [initialPWs;linearCalibrationPWs];
		 nAdd = nAdd - 1;
       else
         linearCalibrationPWs = '';
       end
   catch
       % This will sometimes fail if all of the bagged
       % trials don't run.
       linearCalibrationPWs = '';
   end	   
end

if nAdd > 0
	[nTest, ~] = size(initialPWs);
	nSpreadPerInitial = ceil(nAdd/nTest);
	for addCounter = 1:nTest
		if nAdd > 0
			newPWSets = spreadPWsToNeighbors(initialPWs(addCounter,:), myVPop.coeffsTable, nSpreadPerInitial);
			newPWSets = newPWSets(2:end,:);
			nAddCur = min(nSpreadPerInitial, nAdd);
			nAdd = nAdd - nAddCur;
			initialPWs = [initialPWs;newPWSets(1:nAddCur,:)];
		end
	end
end

if nAddScores > 0
	initialPWs = [initialPWs;vpScores];
end
end