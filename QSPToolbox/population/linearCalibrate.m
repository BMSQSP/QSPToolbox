function [myVPopOptim,optimResults] = linearCalibrate(myVPop,myLinearCalibrateOptions)
% Performs fast linear least-squares optimization of prevalence weights for
% data in a QSPToolbox, using the 'lsqnonneg' function.
%
% ARGUMENTS
%  myVPop 					(VPopRECISTnoBin object)
%  myLinearCalibrateOptions	(linearCalibrateOptions object) Specifies optimization parameters. 
%
% RETURNS
%  myVPopOptim 				(VPop class) VPop with myVPop.pws set to the optimal
% 							 prevalence weights renormalized so that 
%                            they sum to 1, and with the tables
% 							 and goodness-of-fits recalculated.
%  optimResults:			 (struct) Contains information about the optimization
% 							 results. Contains the following fields:
%                            - optimalPrevalenceWeights: (vector) optimal prevalence weights, prior to
%                                                         renormalization
%                            - normOfResidualsSquared:   (scalar) See 'lsqnonneg' documentation
%                            - residuals: 				 (vector): See 'lsqnonneg' documentation
%                            - exitFlag: 				 See 'lsqnonneg' documentation
%                            - output: 					 See 'lsqnonneg' documentation
%                            - lambda: 					 See 'lsqnonneg' documentation
%  fractionOfNonZeroPrevalenceWeights:	(scalar) Fraction of the optimal
% 										prevalence weights that are greater than zero.
% sumOfOptimalPrevalenceWeights: 		(scalar) Sum of the optimal prevalence weights
% 										prior to renormalization
%

% to do:
% - add 2D distributions
% - check whether using sparse computations speeds or slows things down
% - check whether changing tolerance of optimization changes the effN
% - test whether preallocating makes a difference in performance
% - develop iterative prevalence weighting scheme. could use fmincon in an outer loop to
%   optimize the prevalence weight
% additional ideas
% - methods for sampling and alternate solutions, if rank deficient
% - evaluate impact of CDF discretization on spreads the PWs

continueFlag = false;
if nargin > 2
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: myVPop and optionally myLinearCalibrateOptions.'])
    continueFlag = false;
elseif nargin > 1
    continueFlag = true;
elseif nargin > 0
    continueFlag = true;	
	myLinearCalibrateOptions = linearCalibrateOptions;
else
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: myVPop and optionally myLinearCalibrateOptions.'])
end

if continueFlag
    if ~ismember(class(myVPop),{'VPopRECISTnoBin'})
        continueFlag = false;
        warning(['Input VPop not recognized in call to ',mfilename,'.  A VPopRECISTnoBin is needed.'])
    end       
end

if continueFlag

	tictocMain = tic();
	nVPs = size(myVPop.simData.Data,2);

	% The prevalence weights are required to renormalize VP prevalence weights
	% in cases where some VP simulation values are NaN since they dropped off
	% therapy. With VPs missing, the prevalence weights need to be renormalized
	% before calculating the weighted mean and standard deviation of a
	% simulation variable across the virtual population.
	ignoreNaN = false;
	if strcmp(myLinearCalibrateOptions.priorPrevalenceWeightAssumption,'uniform')    
		pws = (1/nVPs)*ones(1,nVPs);
	elseif strcmp(myLinearCalibrateOptions.priorPrevalenceWeightAssumption,'specified')
		pws = myVPop.pws;
	elseif strcmp(myLinearCalibrateOptions.priorPrevalenceWeightAssumption,'ignoreDropout')
		pws = NaN;
		ignoreNaN = true;
	end

	simDataRowInfoTbl = cell2table(myVPop.simData.rowInfo,'VariableNames',myVPop.simData.rowInfoNames);

	% Initialize 'A' in Ax=b
	% Will eventually be matrix of size 'nObservations' x 'nVPs'
	independentVarVal = [];

	% Initialize 'b' in Ax=b
	% Will eventually be vector of size 'nObservations' x 1
	responseVal = [];

	% Each observation is given a weight equal to the 'weight' column in the
	% myVPop data tables.
	% Will eventually be vector of size 'nObservations' x 1
	weightsOfObservations = [];

	% The user has the ability to specify different weights for different data
	% groups (e.g., binTable vs distTable vs etc)
	% Will eventually be vector of size 'nObservations' x 1; thus, this vector
	% will contain redundant information, since different observations of the
	% same group will contain the same value for 'weightsOfDataGroups'
	weightsOfDataGroups = [];

	% String descriptions of the different observations -- useful for
	% debugging.
	% Will eventually be vector of size 'nObservations' x 1
	obsDescriptions = {};

	% Extract data for 'binTable'
	tictocDataType = tic();
	if myLinearCalibrateOptions.verbose, disp('Extracting data related to binTable...'); end
	for iBinTableRow = 1:height(myVPop.binTable)
	% Loop through rows in table
		
		% Determine the index of the row in the simData table that corresponds
		% to this particular 'binTable' row:
		iMaskSimDataRow = ...
			simDataRowInfoTbl.time == myVPop.binTable.time(iBinTableRow) & ...
			strcmp(simDataRowInfoTbl.interventionID,myVPop.binTable.interventionID{iBinTableRow}) & ...
			strcmp(simDataRowInfoTbl.elementID,myVPop.binTable.elementID{iBinTableRow}) & ...
			strcmp(simDataRowInfoTbl.elementType,myVPop.binTable.elementType{iBinTableRow}) & ...
			strcmp(simDataRowInfoTbl.expVarID,myVPop.binTable.expVarID{iBinTableRow});
			
		% extract bin edges
		binTableVariableNames = myVPop.binTable.Properties.VariableNames;
		binEdges = -Inf;
		iBinEdge = 1;
		while true
			variableNameParticular = ['binEdge' num2str(iBinEdge)];
			variableIndex = strcmp(binTableVariableNames,variableNameParticular);
			if ~any(variableIndex)
				break;
			end
			binEdges = [binEdges myVPop.binTable.(variableNameParticular)(iBinTableRow)];
			iBinEdge = iBinEdge + 1;
		end
		binEdges = [binEdges Inf];
		
		% Extract the simulation values, and bin them according to the binEdges
		% determined above
		simVals = myVPop.simData.Data(iMaskSimDataRow,:);
		simValsBinned = arrayfun(@(x) find(binEdges<x,1,'last'),simVals);
		
		% Loop through bins and define an observation row for each bin
		for iBin = 1:length(binEdges)-1
			% The observed response is taken to be the experimentally observed
			% probability of this bin
			responseValParticular = myVPop.binTable.(['expBin' num2str(iBin)])(iBinTableRow);
			% Description of this observation:
			descriptionParticular = ['binTable; Row ' num2str(iBinTableRow) '; Bin ' num2str(iBin)];
			if responseValParticular ~= 0
				% The values for the independent variables will be set to 1 for
				% VPs in the bin, and 0 for those outside the binL
				independentVarValParticular = simValsBinned;
				independentVarValParticular(independentVarValParticular~=iBin) = 0;
				independentVarValParticular(independentVarValParticular==iBin) = 1;
				independentVarVal = [independentVarVal; independentVarValParticular];
				responseVal = [responseVal; responseValParticular];
				weightParticular = myVPop.binTable.weight(iBinTableRow);
				if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
				% normalize so that fit doesn't depend on magnitude
					weightParticular = weightParticular / abs(responseValParticular); 
				end
				weightsOfObservations = [weightsOfObservations weightParticular];
				weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.binTableGroupWeight];
				obsDescriptions = [obsDescriptions; {descriptionParticular}];
			else
				warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
			end
		end
	end
	if myLinearCalibrateOptions.verbose, disp(['Finished extracting data related to binTable. Elapsed time [min]: ' num2str(toc(tictocDataType)/60)]); end

	% Extract data from distTable
	tictocDataType = tic();
	if myLinearCalibrateOptions.verbose, disp('Extracting data related to distTable...'); end
	nRows = height(myVPop.distTable);
	for iDistTableRow = 1:nRows
	% loop through table rows
		tictocRow = tic();
		if myLinearCalibrateOptions.verbose, disp(['Extracting data from distTable Row ' num2str(iDistTableRow) ' of ' num2str(nRows) '...']); end
		% Some of the code below is adapted from the function
		% 'plotDistCDFVPop.m':
		myTable = myVPop.distTable;
		expN = myTable{iDistTableRow,'expN'};
	%     predN = myTable{iDistTableRow,'predN'};
		expW = 1./(expN*ones(1,expN));
		predW = myTable{iDistTableRow,'predProbs'}{1};
	%     expSample = myTable{iDistTableRow,'expSample'}{1};
	%     predSample = myTable{iDistTableRow,'predSample'}{1};
		SC = myTable{iDistTableRow,'combinedPoints'}{1};
		expInd = myTable{iDistTableRow,'expCombinedIndices'}{1};
		predInd = myTable{iDistTableRow,'simCombinedIndices'}{1};		
		[CDFexp, ~] = alignCDFsPreGrid(SC, expInd, predInd, expW, predW);
		% If the user hasn't specified to fit only certain probabilities on the
		% CDF curve, then fit all of the probabilities. This may take a very
		% long time if the CDF contains many points.
		if isempty(myLinearCalibrateOptions.cdfProbsToFit)
			myLinearCalibrateOptions.cdfProbsToFit = CDFexp;
        end
        %cdfVarVal = nan(,nVPs)
		for iIncludedProbs = 1:length(myLinearCalibrateOptions.cdfProbsToFit)
		% loop through the probabilities to fit
			targetProb = myLinearCalibrateOptions.cdfProbsToFit(iIncludedProbs);
			% find the index that corresponds to the closest probability to the
			% target probability:
			closestIndicesToTargetProb = find(abs(CDFexp-targetProb) == min(abs(CDFexp-targetProb)));
			% In case there are multiple closest indices, take the middle of
			% the indices:
			middleClosestIndexToTargetProb = round(closestIndicesToTargetProb(1) + (closestIndicesToTargetProb(end)-closestIndicesToTargetProb(1))/2);
			% The response value is the experimentally observed CDF
			% probability:
			responseValParticular = CDFexp(middleClosestIndexToTargetProb);
			descriptionParticular = ['distTable; Row ' num2str(iDistTableRow) '; Included Prob ' num2str(iIncludedProbs)];
			if responseValParticular ~= 0
				% The values of the independent variables will be set to 1 for
				% VPs that contribute to the particular point on the CDF, and
				% to 0 for VPs that don't contribute to that point:
				if isnan(myVPop.distTable.simCombinedIndices{iDistTableRow}(middleClosestIndexToTargetProb))
					vpIndicesInBin = [];
				else
					vpIndicesInBin = myVPop.distTable.predIndices{iDistTableRow}(1:myVPop.distTable.simCombinedIndices{iDistTableRow}(middleClosestIndexToTargetProb));
				end
				independentVarValParticular = zeros(1,nVPs);
				independentVarValParticular(vpIndicesInBin) = 1;
				independentVarVal = [independentVarVal; independentVarValParticular];
				responseVal = [responseVal; responseValParticular];
				weightParticular = myVPop.distTable.weight(iDistTableRow);
				if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
				% normalize so that fit doesn't depend on magnitude
					weightParticular = weightParticular / abs(responseValParticular); 
				end
				weightsOfObservations = [weightsOfObservations weightParticular];
				weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.distTableGroupWeight];
				obsDescriptions = [obsDescriptions; {descriptionParticular}];
			else
				warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
			end
		end
		tictocRow = tic();
		if myLinearCalibrateOptions.verbose, disp(['Finished extracting data from distTable Row ' num2str(iDistTableRow) '. Elapsed time [min]: ' num2str(toc(tictocRow)/60)]); end
	end
	if myLinearCalibrateOptions.verbose, disp(['Finished extracting data related to distTable. Elapsed time [min]: ' num2str(toc(tictocDataType)/60)]); end

	% Extract data related to the 'brTable':
	% RECIST Class definitions
	% 0 = CR, 1 = PR, 2 = SD, 3 = PD
	recistClassStr = {'CR','PR','SD','PD'};
	resistClassNum = [0 1 2 3];
	tictocDataType = tic();
	if myLinearCalibrateOptions.verbose, disp('Extracting data related to brTableRECIST...'); end
	for iBRTblRow = 1:height(myVPop.brTableRECIST)
	% Loop through table rows
		% Extract the simulation data:
		iMaskSimDataRow = myVPop.simData.brRows == iBRTblRow;
		brDataParticular = myVPop.simData.brData(iMaskSimDataRow,:);
		for iRECISTClass = 1:length(resistClassNum)
		% loop through RECIST classes
			% The response value is the experimentally observed probability for
			% this RECIST class:
			responseValParticular = myVPop.brTableRECIST.(['exp' recistClassStr{iRECISTClass}])(iBRTblRow);
			descriptionParticular = ['brTableRECIST; Row ' num2str(iBRTblRow) '; RECIST Class ' num2str(iRECISTClass)];
			if responseValParticular ~= 0
				% The values of the independent variables will be set to 1 for
				% VPs that are in the RECIST class, and 0 for those that aren't
				vpsInClass = brDataParticular == resistClassNum(iRECISTClass);
				independentVarValParticular = vpsInClass;
				independentVarVal = [independentVarVal; independentVarValParticular];
				responseVal = [responseVal; responseValParticular];
				weightParticular = myVPop.brTableRECIST.weight(iBRTblRow);
				if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
				% normalize so that fit doesn't depend on magnitude
					weightParticular = weightParticular / abs(responseValParticular); 
				end
				weightsOfObservations = [weightsOfObservations weightParticular];
				weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.brTableRECISTGroupWeight];
				obsDescriptions = [obsDescriptions; {descriptionParticular}];
			else
				warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
			end
		end
	end
	if myLinearCalibrateOptions.verbose, disp(['Finished extracting data related to brTableRECIST. Elapsed time [min]: ' num2str(toc(tictocDataType)/60)]); end

	tictocDataType = tic();
	if myLinearCalibrateOptions.verbose, disp('Extracting data related to rTableRECIST...'); end
	for iRTblRow = 1:height(myVPop.rTableRECIST)
		iMaskSimDataRow = find(myVPop.simData.rRows == iRTblRow);
		rDataParticular = myVPop.simData.rData(iMaskSimDataRow,:);
		% 0 = CR, 1 = PR, 2 = SD, 3 = PD
		recistClassStr = {'CR','PR','SD','PD'};
		resistClassNum = [0 1 2 3];
		for iRECISTClass = 1:length(resistClassNum)
		% loop through RECIST classes
			% The response value is the experimentally observed probability for
			% this RECIST class:
			responseValParticular = myVPop.rTableRECIST.(['exp' recistClassStr{iRECISTClass}])(iRTblRow);
			descriptionParticular = ['rTableRECIST; Row ' num2str(iRTblRow) '; RECIST Class ' num2str(iRECISTClass)];
			if responseValParticular ~= 0
				% The values of the independent variables will be set to 1 for
				% VPs that are in the RECIST class, and 0 for those that aren't
				vpsInClass = rDataParticular == resistClassNum(iRECISTClass);
				independentVarValParticular = vpsInClass;
				independentVarVal = [independentVarVal; independentVarValParticular];
				responseVal = [responseVal; responseValParticular];
				weightParticular = myVPop.rTableRECIST.weight(iRTblRow);
				if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
				% normalize so that fit doesn't depend on magnitude
					weightParticular = weightParticular / abs(responseValParticular); 
				end
				weightsOfObservations = [weightsOfObservations weightParticular];
				weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.rTableRECISTGroupWeight];
				obsDescriptions = [obsDescriptions; {descriptionParticular}];
			else
				warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
			end
		end
	end
	if myLinearCalibrateOptions.verbose, disp(['Finished extracting data related to rTableRECIST. Elapsed time [min]: ' num2str(toc(tictocDataType)/60)]); end

	% Extract data for mean and standard deviation
	tictocDataType = tic();
	if myLinearCalibrateOptions.verbose, disp('Extracting data related to mnSDTable...'); end
	for iSumStatRow = 1:height(myVPop.mnSDTable)
	% Loop through table rows
		% Extract the appropriate simulation data:
		iMaskSimDataRow = ...
			simDataRowInfoTbl.time == myVPop.mnSDTable.time(iSumStatRow) & ...
			strcmp(simDataRowInfoTbl.interventionID,myVPop.mnSDTable.interventionID{iSumStatRow}) & ...
			strcmp(simDataRowInfoTbl.elementID,myVPop.mnSDTable.elementID{iSumStatRow}) & ...
			strcmp(simDataRowInfoTbl.elementType,myVPop.mnSDTable.elementType{iSumStatRow}) & ...
			strcmp(simDataRowInfoTbl.expVarID,myVPop.mnSDTable.expVarID{iSumStatRow});
		simVals = myVPop.simData.Data(iMaskSimDataRow,:);
		
		% Some simVals will be NaN for VPs that dropped off of therapy. Thus,
		% we need to scale the prevalence weights so that the sum of the
		% weights for the VPs still on therapy adds to 1.
		nanValsIndexMask = isnan(simVals);
		if ignoreNaN
			if any(nanValsIndexMask)
				continue;
			end
			% If we reached here, there are no NaN VPs and thus sumPWsForNonNanVPs = 1
			sumPWsForNonNanVPs = 1;
		else
			sumPWsForNonNanVPs = sum(pws(~nanValsIndexMask));
		end
		
		% Fit mean
		% The response value is the experimental mean:
		responseValParticular = myVPop.mnSDTable.expMean(iSumStatRow);
		descriptionParticular = ['mnSDTable; Row ' num2str(iSumStatRow) '; mean'];
		if responseValParticular ~= 0
			% The values for the independent variables are simply the
			% simulation values
			independentVarValParticular = simVals;
			% renormalize prevalence weights wtih NaN VPs purged. We will do
			% this by scaling the values for the independent variables:
			independentVarValParticular = (1/sumPWsForNonNanVPs)*independentVarValParticular; 
			% NaN values will cause issues in the optimization algorithm. 
			% Since we renormalized the prevalence weight, change NaN values to
			% zero (which shouldn't affect the calculation):
			independentVarValParticular(nanValsIndexMask) = 0;
			independentVarVal = [independentVarVal; independentVarValParticular];
			responseVal = [responseVal; responseValParticular];
			weightParticular = myVPop.mnSDTable.weightMean(iSumStatRow);
			if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
			% normalize so that fit doesn't depend on magnitude
				weightParticular = weightParticular / abs(responseValParticular); 
			end
			weightsOfObservations = [weightsOfObservations weightParticular];
			weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.mnSDTableGroupWeight];
			obsDescriptions = [obsDescriptions; {descriptionParticular}];
		else
			warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
		end
		
		% Fit SD (actually, variance)
		% The response value will be the variance:
		responseValParticular = myVPop.mnSDTable.expSD(iSumStatRow)^2;
		descriptionParticular = ['mnSDTable; Row ' num2str(iSumStatRow) '; stdev'];
		if responseValParticular ~= 0
			% The values for the independent variables are the individual
			% variances (i.e., 'residualsSquared'). We cannot yet know the
			% predicted mean, since that would require prior knowledge of the
			% optimal prevalence weights. But, since we are fitting the mean,
			% let's assume that the predicted mean will be (approximately)
			% equal to the experimental mean.
			expMean = myVPop.mnSDTable.expMean(iSumStatRow);
			residualsSquared = (simVals-expMean).^2;
			independentVarValParticular = residualsSquared;
			% renormalize prevalence weights wtih NaN VPs purged. We will do
			% this by scaling the values for the independent variables:
			independentVarValParticular = (1/sumPWsForNonNanVPs)*independentVarValParticular; 
			% NaN values will cause issues in the optimization algorithm. 
			% Since we renormalized the prevalence weight, change NaN values to
			% zero (which shouldn't affect the calculation):
			independentVarValParticular(nanValsIndexMask) = 0;
			independentVarVal = [independentVarVal; independentVarValParticular];
			responseVal = [responseVal; responseValParticular];
			weightParticular = myVPop.mnSDTable.weightSD(iSumStatRow);
			if strcmp(myLinearCalibrateOptions.responseValTransformation,'relative')
			% normalize so that fit doesn't depend on magnitude
				weightParticular = weightParticular / abs(responseValParticular); 
			end
			weightsOfObservations = [weightsOfObservations weightParticular];
			weightsOfDataGroups = [weightsOfDataGroups myLinearCalibrateOptions.mnSDTableGroupWeight];
			obsDescriptions = [obsDescriptions; {['mnSDTable; Row ' num2str(iSumStatRow) '; stdev']}];
		else
			warning(['Ignoring responses that are equal to zero, since the least-squares fitting is performed on a relative scale so that it does not depend on absolute values. Thus, the following observation is being ignored: ' descriptionParticular]);
		end
	end
	if myLinearCalibrateOptions.verbose, disp(['Finished extracting data related to mnSDTable. Elapsed time [min]: ' num2str(toc(tictocDataType)/60)]); end

	%     nIndependentVars = size(independentVarVal,2);
	%     
	%     nonNegativityConstraints_A = diag(-1*ones(1,nIndependentVars),0);
	%     nonNegativityConstraints_b = zeros(nIndependentVars,1);
	% [x,resnorm,residual,exitflag,output,lambda] = lsqlin(independentVarVal,responseVal,nonNegativityConstraints_A,nonNegativityConstraints_b);

	% Put variables into pre-allocated form before running fit:
	tictocPreallocate = tic();
	if myLinearCalibrateOptions.verbose, disp('Preallocating matrices...'); end
	independentVarValPreAllocated = independentVarVal;
	clear('independentVarVal');
	responseValPreAllocated = responseVal;
	clear('responseVal');
	if myLinearCalibrateOptions.verbose, disp(['Finished preallocating matrices. Elapsed time [min]: ' num2str(toc(tictocPreallocate)/60)]); end

	% factor in weights
	tictocWeights = tic();
	if myLinearCalibrateOptions.verbose, disp('Multiplying weights...'); end
	% lsqWeights = sparse(diag(sqrt(weightsOfObservations.*weightsOfDataGroups)));
	% lsqWeights = sqrt(weightsOfObservations.*weightsOfDataGroups);
	lsqWeights = weightsOfObservations.*weightsOfDataGroups;
	lsqWeights = spdiags(lsqWeights',0,length(lsqWeights),length(lsqWeights));
	independentVarValPreAllocated = lsqWeights*independentVarValPreAllocated;
	responseValPreAllocated = lsqWeights*responseValPreAllocated;
	if myLinearCalibrateOptions.verbose, disp(['Finished multiplying weights. Elapsed time [min]: ' num2str(toc(tictocWeights)/60)]); end

	% Run the optimization
	tictocOptimization = tic();
	if myLinearCalibrateOptions.verbose, disp('Performing optimization...'); end
	lsqnonnegOptions = optimset;

    % Default tolerance from 'lsqnonneg' documentation:
    defaultTol = 10*eps*norm(independentVarValPreAllocated,1)*length(independentVarValPreAllocated);
	tol = defaultTol*myLinearCalibrateOptions.toleranceScalingFactor;
	lsqnonnegOptions.TolX = tol;
	
	[optimResults.optimalPrevalenceWeights,optimResults.normOfresidualsSquared,optimResults.residuals,optimResults.exitFlag,optimResults.output,optimResults.lambda] = lsqnonneg(independentVarValPreAllocated,responseValPreAllocated,lsqnonnegOptions);
	if myLinearCalibrateOptions.verbose, disp(['Finished performing optimization. Elapsed time [min]: ' num2str(toc(tictocOptimization)/60)]); end

	if myLinearCalibrateOptions.verbose, disp(['Finished program. Total elapsed time [min]: ' num2str(toc(tictocMain)/60)]); end

	optimResults.sumOfOptimalPrevalenceWeights = sum(optimResults.optimalPrevalenceWeights);

	% Update VPop
	newPWs = optimResults.optimalPrevalenceWeights./sum(optimResults.optimalPrevalenceWeights);

	% Sanity check:
	if abs(sum(newPWs)-1) > 1e-12
		error('Prevalence weights do not sum to 1 after renormalization.')
	end

	myVPopOptim = myVPop;
	myVPopOptim.pws = newPWs';

	% recalculate the tables and goodness of fit 
	myVPopOptim=addPredTableVals(myVPopOptim);
	myVPopOptim=evaluateGOF(myVPopOptim);

	if myLinearCalibrateOptions.createPlotsOfInputVPop
		plotBinVPop(myVPop);
		sgtitle('Input VPop');
		plotBRVPop(myVPop);
		sgtitle('Input VPop');
		plotMnSDVPop(myVPop);
		sgtitle('Input VPop');
		plotDistCDFVPop(myVPop);
		sgtitle('Input VPop');
	end

	if myLinearCalibrateOptions.createPlotsOfOptimizedVPop
		plotBinVPop(myVPopOptim);
		sgtitle('Optimized VPop');
		plotBRVPop(myVPopOptim);
		sgtitle('Optimized VPop');
		plotMnSDVPop(myVPopOptim);
		sgtitle('Optimized VPop');
		plotDistCDFVPop(myVPopOptim);
		sgtitle('Optimized VPop');
	end

	% scatter plot of new vs old prevalence weights
	if myLinearCalibrateOptions.createPlotsOfInputVPop && myLinearCalibrateOptions.createPlotsOfOptimizedVPop
		figure;
		scatter(myVPop.pws,myVPopOptim.pws);
		xlabel('prevalence weights in input VPop');
		ylabel('optimal prevalence weights');
		title('prevalence weights of optimized vs input VPops');
	end

	% fraction of prevalence weights that are greater than zero
	optimResults.fractionOfNonZeroPrevalenceWeights = sum(myVPopOptim.pws>0)/length(myVPopOptim.pws);

else
    myVPopOptim = VPop;
	optimResults = struct;
    warning(['Unable to run ',mfilename,'.  Returning input VPop and empty optimResults structure.'])
end