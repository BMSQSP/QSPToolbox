classdef LinearCalibration
% Performs fast linear least-squares optimization of prevalence weights for
% data in a QSPToolbox, using the 'lsqlin' or 'lsqnonneg' function. Can
% also perform bootstrapping to calculate confidence intervals, and bagging
% to reduce variability in the calibration.
%
% INPUT ARGUMENTS
% -myVPop: (VPop object) Doesn't need to include the prevalence weights.
%
% PROPERTIES
% -OptimOptions: (LinearCalibrationObtions object) Specifies optimization parameters. 
% Refer to the documentation in LinearCalibrationObtions.m
% -InputVPop: (VPop object) Equal to the 'myVPop' supplied by the user
% -LinearProblemMatrices: (struct) Contains matrices and information
% related to the weighted linear least-squares fitting problem w^(1/2)Ax=w^(1/2)b, where w is
% weights, A is the values of the independent variables, x is the
% prevalence weights we are optimizing, and b is the observation values.
% Contains the following fields:
% --independentVarVals: (matrix) Corresponds to the matrix 'A' in the linear
% problem. nObservations x nVPs
% --observationVals: (vector) Corresponds to the vector 'b' in the linear
% problem. nObservations x 1
% --observationWeights: (vector) Corresponds to the weights 'w' in the linear
% problem. nObservations x 1
% --dataGroup: (vector) The data group index of each observation (a "data
% group" is a set of non-independent observations). nObservations x 1
% --observationDescriptions: (cell array) String descriptions of the 
% different observations -- useful for debugging. nObservations x 1
% --independentVarValsWeighted: (matrix) Corresponds to the product 'w^(1/2)A' in the linear
% problem. nObservations x nVPs
% --observationValsWeighted: (vector) Corresponds to the product 'w^(1/2)b' in the linear
% problem. nObservations x 1
% -OptimalPrevalenceWeightsNormalized: (vector) Prevalence weights returned
% by the optimization, normalized so they sum to 1.
% -OptimalPrevalenceWeightsNormalizedConfidenceIntervals (matrix): 95%
% confidence intervals on 'OptimalPrevalenceWeightsNormalized' 
% -OptimizedVPop: (VPop object) 'InputVPop' with 'InputVPop.pws' set to the optimal
% prevalence weights renormalized so that they sum to 1, and with the tables
% and goodness-of-fits recalculated.
% -IgnoredObservationsDescriptions: (cell array) Descriptions of
% observations which were ignored in the optimization, typically because
% the observation is equal to zero and the fit is performed on a relative
% basis, thus the observation must be ignored so as not to divide by zero.
% -InputVPopStatistics: (struct) Residuals and root-mean-squared error
% based on the prevalence weights of the input VPop.
% -RMSE: (scalar) Root-mean-squared error from the fitting (i.e., on the
% training data set)
% -FractionRunsConverged: (scalar) Fraction of runs which converged
% !!! should make first character lower case
properties
    OptimOptions
    InputVPop
    LinearProblemMatrices
    OptimalPrevalenceWeightsNormalized
    OptimalPrevalenceWeightsNormalizedConfidenceIntervals
    OptimizedVPop
    IgnoredObservationsDescriptions
    InputVPopStatistics
    RMSE
    FractionRunsConverged
end

% Hidden properties:
properties (Hidden = true)
    NumVPs % number of VPs
    TimeElapsedMinutes
    LinearProblemMatricesParticular % Modified LinearProblemMatrices in cases where we want to select particular data groups, etc. This is the matrices which is actually being fitted.
    OptimizationResults % Optimization results
    CrossValidationResults
end

% Private properties:
properties (Access = private)

end


%% PUBLIC METHODS

methods
    function Obj = LinearCalibration(myVPop,varargin)
    % Constructor for initializing a LinearCalibration object.
    %
    % INPUT ARGUMENTS
    % -myVPop: (VPop object) Used for extracting experimental and
    % simulation data. The prevalence weights do not need to be
    % specified.
    % -optimOptions: (LinearCalibrationOptions object) Optional argument
    % for customizing the optimization. Specified as a name-value pair.
    % 
    
        % Parse the optional function arguments:
        p = inputParser;
        addOptional(p,'optimOptions',LinearCalibrationOptions)
        parse(p,varargin{:});
        Obj.OptimOptions = p.Results.optimOptions;
        
        % Save the user-supplied VPop:
        Obj.InputVPop = myVPop;
        
        % Number of VPs:
        Obj.NumVPs = size(Obj.InputVPop.simData.Data,2);
        
        % add path to nnls function
        % Change the current folder to the folder of this m-file.
        % https://www.mathworks.com/matlabcentral/answers/72309-to-change-current-folder-to-the-folder-of-m-file#answer_82461
        dirMem = pwd();
        if(~isdeployed)
          cd(fileparts(which(mfilename)));
        end
        nnlsDir = 'nnls';
        exportFigDir = 'altmany-export_fig-b1a7288';
        if exist(nnlsDir,'dir')
            addpath(genpath(nnlsDir));
        end
        if exist(exportFigDir,'dir')
            addpath(genpath(exportFigDir));
        end
        if(~isdeployed)
        	cd(dirMem);
        end
    end
    
    function Obj = run(Obj,varargin)
    % Constructs the linear matrices, runs optimization, and constructs a
    % new VPop post-optimization.
    %

        timerTotal = tic();

        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'closeParallelPoolWhenFinished',true,@islogical);
        parse(p,varargin{:});
        
        % reset the optimization results
        Obj.OptimizationResults = [];
        
        % Construct matrices:
        if isempty(Obj.LinearProblemMatrices)
            timerMatrixConstruction = tic();
            Obj = Obj.constructLinearProblemMatrices();
            Obj.TimeElapsedMinutes.matrixConstruction = toc(timerMatrixConstruction)/60;
        end
        
        % Run optimization:
        % Start timer:
        timerOptimization = tic();
        % Array of data group indices which we can sample from:
        dataGroups = (unique(Obj.LinearProblemMatrices.dataGroup))';
        % Number of data groups:
        nDataGroups = length(dataGroups);
        % Run optimization based on the specified method:
        if strcmpi(Obj.OptimOptions.method,"bestFit")
        % Run a single best-fit optimization
            Obj = Obj.runOptimization();
            Obj.FractionRunsConverged = 0;
            if ~any(isnan(Obj.OptimizationResults.optimalPrevalenceWeightsNormalized))
                Obj.FractionRunsConverged = 1;
            end
                
        elseif strcmpi(Obj.OptimOptions.method,"bootstrap")
        % Run a single best-fit optimization, and then calculate the
        % confidence intervals on the prevalence weights using
        % bootstrapping
            % First, run the single best-fit optimization:
            Obj = Obj.runOptimization();
            
            % Specify the options for bootstrapping:
            bootstrapOptions = statset();
            bootstrapOptions.UseParallel = true;
            bootstapFunctionHandle = @(includedDataGroupsParticular) bootstrapFunction(Obj,includedDataGroupsParticular);
            
            % Calculate the confidence intervals using bootstrapping:
            Obj.OptimalPrevalenceWeightsNormalizedConfidenceIntervals = bootci(Obj.OptimOptions.nBootstrapIterations,{bootstapFunctionHandle,dataGroups},'Options',bootstrapOptions);

        elseif strcmpi(Obj.OptimOptions.method,"bagging")
        % Use bagging to reduce the variability in the estimation of the
        % optimal prevalence weights
            
            % Create a vector of VP indices that 
            vpIndices = (1:Obj.NumVPs)';

            % Convert the fraction of VPs per iteration into the number of VPs
            nVPsPerBaggingIteration = floor(Obj.OptimOptions.fractionVPsPerBaggingIteration*Obj.NumVPs);

            % Perform the bagging
            % First, initialize a cell array to store the VP indices
            % included in each iteration:
            includedVPIndices = cell(1,Obj.OptimOptions.nBootstrapIterations);
            % Also initialize a cell array to store the optimal
            % prevalence weights calculated in each iteration:
            optimizationPrevalenceWeightsOrderedAccordingToSampledVPIndices = cell(1,Obj.OptimOptions.nBootstrapIterations);
            runConverged = zeros(1,Obj.OptimOptions.nBootstrapIterations);
            parfor iBagging = 1:Obj.OptimOptions.nBootstrapIterations % parfor
                % In case we want to be verbose:
                % tictoc = tic();
                % disp(['** BAGGING: starting Iteration ' num2str(iBootstrap) ' of ' num2str(Obj.OptimOptions.nBootstrapIterations) ' **']);
                
                % Sample the data groups for this iteration. Like bootstrapping, the sampling is
                % performed with replacement, so that some data groups may
                % be represented more than once (in which case they would
                % be weighted more highly in the fit).
                includedDataGroups = datasample(dataGroups,nDataGroups);
                
                % Sample the VPs for this iteration, without replacement
                % (since having a VP represented more than once can create
                % issues in the fit).
                includedVPIndices{iBagging} = datasample(vpIndices,nVPsPerBaggingIteration,'Replace',false);
                
                % Rebuild the matrices based on the sampled VPs and data
                % groups:
                ObjParticular = Obj.selectFromLinearProblemMatrices('includedDataGroups',includedDataGroups,'includedVPIndices',includedVPIndices{iBagging});
                
                % Run the optimization:
                ObjParticular = ObjParticular.runOptimization();
                
                % Save the optimal prevalence weights for this iteration.
                % Note that the results for each iteration is ordered
                % according to the ordering of the sampled VP indices
                % (hence why above we also saved the ordering of the VP indices)
                optimizationPrevalenceWeightsOrderedAccordingToSampledVPIndices{iBagging} = ObjParticular.OptimizationResults.optimalPrevalenceWeightsNormalized;
                
                if ~any(isnan(ObjParticular.OptimizationResults.optimalPrevalenceWeightsNormalized))
                    runConverged(iBagging) = 1;
                end
                
                % In case we want to be verbose:
                % disp(['Time elapsed for Bagging Iteration ' num2str(iBootstrap) ' [seconds]: ' num2str(toc(tictoc))]);
            end

            Obj.FractionRunsConverged = sum(runConverged)/Obj.OptimOptions.nBootstrapIterations;
            
            % Consolidate the optimal prevalence weights into a matrix.
            % The matrix is initialized with zeros so that any VP not
            % included in an iteration is automatically assigned a
            % prevalence weight of zero for that iteration (this is what 
            % we want, so we'll keep it that way). 
            optimalPrevalenceWeightsByIteration = zeros(Obj.NumVPs,Obj.OptimOptions.nBootstrapIterations);
            for iBagging = 1:Obj.OptimOptions.nBootstrapIterations
                % The indexing of the left-hand side of the following line
                % is performed to account for the fact that the ordering of
                % the VPs is different in different iterations of the
                % bagging:
                optimalPrevalenceWeightsByIteration(includedVPIndices{iBagging},iBagging) = optimizationPrevalenceWeightsOrderedAccordingToSampledVPIndices{iBagging};
            end

            % Theory suggests that the mean of the predictions of each
            % iteration is the same as the prediction from the mean
            % prevalence weight. Thus, we can save only the mean
            % prevalence weight:
            Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization = nanmean(optimalPrevalenceWeightsByIteration,2);
        
            % Normalize the optimal prevalence weights:
            Obj.OptimizationResults.optimalPrevalenceWeightsNormalized = ...
                Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization./sum(Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization);
        else
            error('Unrecognized optimization method.');
        end
        % Post-optimization analysis:
        Obj.OptimalPrevalenceWeightsNormalized = Obj.OptimizationResults.optimalPrevalenceWeightsNormalized;
        Obj = Obj.performPostOptimizationAnalysis();
        Obj.RMSE = Obj.OptimizationResults.rmseConsideringOptimalPrevalenceWeightsNormalized;
        
        % delete the parallel pool unless specified not to do so
        if p.Results.closeParallelPoolWhenFinished
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
        
        Obj.TimeElapsedMinutes.optimization = toc(timerOptimization)/60;
        
        % Construct new VPop according to the optimal prevalence weights
        timerVPopUpdate = tic();
        Obj = Obj.constructNewVPopPostOptimization();
        Obj.TimeElapsedMinutes.vPopUpdate = toc(timerVPopUpdate)/60;
        
        Obj.TimeElapsedMinutes.totalTime = toc(timerTotal)/60;
        
        % Nested function used in the 'bootstrap' method:
        function bestFitPrevalenceWeights = bootstrapFunction(Obj,includedDataGroupsParticular)
            % Rebuild the matrices based on the sampled data groups:
            ObjParticular = Obj.selectFromLinearProblemMatrices('includedDataGroups',includedDataGroupsParticular);
            % Run optimization:
            ObjParticular = ObjParticular.runOptimization();
            % Output the best-fit prevalence weights:
            bestFitPrevalenceWeights = ObjParticular.OptimizationResults.optimalPrevalenceWeightsNormalized;
        end
    end
    
    function Obj = runCrossValidation(Obj)
    % Runs leave-one-out cross-validation. Each iteration of the
    % cross-validation leaves out the observations comprising one data
    % group, fits the remainder of the observations, and tests the fit
    % against the observations in the data group which was left out. Thus,
    % there will be as many iterations as there are data groups.
    %
    
        % Array of data group indices which we can sample from:
        dataGroupIndices = unique(Obj.LinearProblemMatrices.dataGroup)';
        
        % Number of data groups:
        nDataGroups = length(dataGroupIndices);
        
        % The total number of observations:
        nObservations = length(Obj.LinearProblemMatrices.observationVals);
        
        % Initialize a matrix that will store the training residuals from all
        % 'nDataGroups' iterations:
        trainingResidualsSquared = nan(nObservations,nDataGroups);
        
        % Initialize a vector that will store the test residual for each
        % observation (each observation will be tested in exactly one
        % iteration -- when the corresponding data group is the one being
        % left out).
        testResidualsSquared = nan(nObservations,1);
        
        % Loop through iterations and perform cross-validation.
        for iDataGroup = 1:nDataGroups
            dataGroup = dataGroupIndices(iDataGroup);
            tictoc = tic();
            disp(['CROSS-VALIDATION: Testing Data Group Index ' num2str(iDataGroup) ' of ' num2str(nDataGroups)]);
            
            % Masking index for selecting observations in the data group:
            iMaskDataGroup = Obj.LinearProblemMatrices.dataGroup == dataGroup;
            
            % The data groups to include in the fit for this iteration:
            dataGroupsToInclude = dataGroupIndices;
            dataGroupsToInclude(iDataGroup) = [];
            
            % Isolate the observations corresponding to
            % 'dataGroupsToInclude':
            ObjParticular = Obj.selectFromLinearProblemMatrices('includedDataGroups',dataGroupsToInclude);
            ObjParticular.LinearProblemMatrices = ObjParticular.LinearProblemMatricesParticular;
            
            % Run the fit:
            ObjParticular = ObjParticular.run('closeParallelPoolWhenFinished',false);
            
            % Calculate the training residuals:
            [trainingResiduals,~,~] = ObjParticular.calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs();
            trainingResidualsSquared(~iMaskDataGroup,iDataGroup) = trainingResiduals.^2;
            
            % Calculate the test residuals:
            ObjParticular.LinearProblemMatrices = Obj.LinearProblemMatrices;
            [testResiduals,~,~] = ObjParticular.evaluateDataGroupResponse(dataGroup);
            testResidualsSquared(iMaskDataGroup) = testResiduals.^2;
            
            disp(['Time elapsed for Cross-Validating Data Group Index ' num2str(iDataGroup) ' [seconds]: ' num2str(toc(tictoc))]);
        end
        % Delete fields which don't make sense for cross-validation:
        Obj.OptimizedVPop = [];
        Obj.OptimalPrevalenceWeightsNormalized = [];
        Obj.RMSE = [];
        % Calculate the training and test RMSEs, using weighted means (note that the mean is weighted since the squared residuals are weighted with weights which sum to 1):
        % For the training residuals, we should weight each of them by
        % nDataGroups-1 since this is how many times they are represented.
        % NaN values are expected where the observation was left out in the
        % test set.
        trainingMSE = nansum(trainingResidualsSquared(:)/(nDataGroups-1));
        testMSE = sum(testResidualsSquared);
        Obj.CrossValidationResults.trainingRMSE = sqrt(trainingMSE);
        Obj.CrossValidationResults.testRMSE = sqrt(testMSE);
        
        % close the parallel pool
        poolobj = gcp('nocreate');
        delete(poolobj);
        
        % for debugging:
        % detect which data groups have smaller test residuals compared to
        % training residuals
%         trainingResidualsSquaredWeighted = trainingWeights.*trainingResidualsSquared./sum(trainingWeights(:));
%         testResidualsSquaredWeighted = testWeights.*testResidualsSquared./sum(testWeights);
%         trainingResidualsSquaredWeightedConsolidatedPerObservation = sum(trainingResidualsSquaredWeighted,2);
%         testMinusTrainingResidualsSquaredWeighted = testResidualsSquaredWeighted-trainingResidualsSquaredWeighted;
%         crossValidationResidualsTable = table((1:nObservations)',Obj.LinearProblemMatrices.observationDescriptions,Obj.LinearProblemMatrices.dataGroup,trainingResidualsSquaredWeightedConsolidatedPerObservation,testResidualsSquaredWeighted,testMinusTrainingResidualsSquaredWeighted,...
%             'VariableNames',{'observationIndex','observationDescriptions','dataGroup','trainingResidualsSquaredWeighted','testResidualsSquaredWeighted','testMinusTrainingResidualsSquaredWeighted'});
%         crossValidationResidualsSortedTable = sortrows(crossValidationResidualsTable,'testMinusTrainingResidualsSquaredWeighted');
    end  

    function resultsOut = useCrossValidationToOptimizeBagging(Obj,outputDir,varargin)
    % Guides the optimization of the fraction of VPs to be selected per
    % bagging iteration. This function will generate plots of
    % root-mean-squared error (RMSE), effective N, and run time; as a
    % function of the fraction, and will save these plots in 'outputDir'.
    %
    % INPUT ARGUMENTS
    % -outputDir: (character array) Directory where to store the results.
    % Must end with '/'.
    % -vPopToPredictAgainst: (VPop object) If an external prediction is
    % desired, supply the data to predict against using this argument.
    % Specified as a name-value pair. Default is [], in which case no
    % external prediction is performed.
    % -fractionVPsPerBaggingIterationToAnalyze: (numeric vector)Indicates
    % which fractions to analyze. This argument is optional, and is
    % specified as a name-value pair. The default is '0.05:0.05:1'.
    % -calculateGOF: (scalar boolean) Indicates whether the function should
    % calculare and return the composite goodness-of-fit calculated using
    % the QSP Toolbox. This argument is optional, and is
    % specified as a name-value pair. The default is 'false'.
    %
    % OUTPUTS
    % -resultsOut: (struct) Contains the following fields:
    % --optimizationResultsTable (table): Contains the cross-validation
    % results for each of the fractions
    % --inputVPop: (struct) Contains results from analysis of
    % 'Obj.InputVPop', in case we'd like to use them as a comparison in
    % order to judge the optimization
    % -resultsTable: (table) Results of analysis.
    %
        
        tictocMain = tic();
        
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'vPopToPredictAgainst',[]);
        addOptional(p,'fractionVPsPerBaggingIterationToAnalyze',[0.001 0.005 0.01 0.05 0.1 0.25 0.5 0.75 1],@isnumeric);
        % addOptional(p,'fractionVPsPerBaggingIterationToAnalyze',[0.01:0.01:0.04, 0.05:0.05:1],@isnumeric);
        addOptional(p,'calculateGOF',false);
        parse(p,varargin{:});
        
        % Create the output directory if it doesn't already exist
        if ~exist(outputDir,'dir')
            mkdir(outputDir);
        end

        % Specify to use bagging:
        Obj.OptimOptions.method = "bagging";
        
        % Start a log file for this run:
        fid = fopen([outputDir 'runLog.txt'],'w');

        % The number of fractions to test:
        nFractions = length(p.Results.fractionVPsPerBaggingIterationToAnalyze);
        
        % Let's do a simple best-fit (with no bagging) run, so that we can
        % compare against it:
        ObjNoBagging = Obj;
        ObjNoBagging.OptimOptions.method = "bestFit";
        outputDirParticular = [outputDir 'noBagging/'];
        [results.noBaggingRMSE,results.noBaggingEffN,results.noBaggingGOF] = runCrossValidationWithComparisons(ObjNoBagging,outputDirParticular,p);
        
        % If adequate information is speified in the InputVPop, calculate
        % some statistics for comparison:
        if ~isempty(Obj.InputVPop.pws)
            % Effective N:
            results.inputVPopEffN.effNWholeFit = Obj.calculateEffN(Obj.InputVPop.pws);
            % RMSE:
            Obj = Obj.calculateInputVPopStatistics('vPopToPredictAgainst',p.Results.vPopToPredictAgainst,'saveFigDir',[outputDir '/inputVPopPredictions/'],'closeFigsAutomatically',true);
            results.inputVPopRMSE.rmseWholeFit = Obj.InputVPopStatistics.rootMeanSquaredError;
            if ~isempty(p.Results.vPopToPredictAgainst)
                results.inputVPopRMSE.rmsePrediction = Obj.InputVPopStatistics.predictionRootMeanSquaredError;
            end
            % Goodness of fit:
            if p.Results.calculateGOF && ~isempty(Obj.InputVPop.gof)
                results.inputVPopRMSE.gofWholeFit = Obj.InputVPop.gof;
            end
        else
            results.inputVPop = [];
        end
        
        % Loop through the specified bagging fractions, and perform
        % cross-validation for each of the freactions:
        for iFraction = 1:nFractions
            tictocParticularFraction = tic();
            
            disp(['------------- STARTING CROSS-VALIDATION FOR FRACTION INDEX ' num2str(iFraction) ' of ' num2str(nFractions) ' -------------']);
            
            % Specify the fraction VPs parameter for this iteration:
            Obj.OptimOptions.fractionVPsPerBaggingIteration = p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction);
            
            % Make a directory specific to this fraction:
            fractionVPsStr = strrep(num2str(p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction)),'.','p');
            outputDirParticular = [outputDir 'fractionIndex_' num2str(iFraction) '__fraction_' fractionVPsStr '/'];

            % run the cross-validation:
            [results.baggingRMSE(iFraction),results.baggingEffN(iFraction),results.baggingGOF(iFraction),results.baggingRunTime(iFraction)] = runCrossValidationWithComparisons(Obj,outputDirParticular,p);
            

            
            disp(['Time elapsed [minutes] for Cross-Validation Fraction Index ' num2str(iFraction) ': ' num2str(toc(tictocParticularFraction)/60)]);
        end
        for iFraction = 1:nFractions
            results.baggingRMSE(iFraction).fractionVPsPerBaggingIteration = p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction);
            results.baggingEffN(iFraction).fractionVPsPerBaggingIteration = p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction);
            results.baggingGOF(iFraction).fractionVPsPerBaggingIteration = p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction);
            results.baggingRunTime(iFraction).fractionVPsPerBaggingIteration = p.Results.fractionVPsPerBaggingIterationToAnalyze(iFraction);
        end

        % Close the parallel pool
        poolobj = gcp('nocreate');
        delete(poolobj);
        
        % Put the results into a table
        resultsFieldNamesTable = fields(results);
        for iFieldTable = 1:length(resultsFieldNamesTable)
            resultsFieldNamesCols = fields(results.(resultsFieldNamesTable{iFieldTable}));
            resultsOut.(resultsFieldNamesTable{iFieldTable}) = table();
            for iFieldCols = 1:length(resultsFieldNamesCols)
                resultsTableParticular = ...
                    table([results.(resultsFieldNamesTable{iFieldTable}).(resultsFieldNamesCols{iFieldCols})]','VariableNames',resultsFieldNamesCols(iFieldCols));
                resultsOut.(resultsFieldNamesTable{iFieldTable}) = [resultsOut.(resultsFieldNamesTable{iFieldTable}) resultsTableParticular];
            end
        end

        % Save the results to the output directory
        save([outputDir 'results.mat'],'resultsOut');
        
        % Generate the plots
        generateAndSaveCrossValidationPlots(outputDir,resultsOut,p);
        
        totalTimeElapsedHours = toc(tictocMain)/3600;
        disp(['Total time [hours] elapsed for the bagging optimization: ' num2str(totalTimeElapsedHours)]);
        fprintf(fid,'Total time [hours] elapsed for the bagging optimization: %f.\n\n', totalTimeElapsedHours);

        % internal function for running cross-validation:
        function [resultsRMSE,resultsEffN,resultsGOF,resultsRunTime] = runCrossValidationWithComparisons(Obj,outputDir,p)
            
            if ~exist(outputDir,'dir')
                mkdir(outputDir);
            end
            
            if isempty(gcp('nocreate'))
                parpool();
            end
            
            % First, run a fit on all of the data
            tictocRun = tic();
            Obj = Obj.run('closeParallelPoolWhenFinished',false);
            resultsRunTime.calibrationTimeForSingleNonCVRunSeconds = toc(tictocRun);
            resultsRMSE.rmseWholeFit = Obj.RMSE;
            resultsEffN.effNWholeFit = Obj.calculateEffN(Obj.OptimizedVPop.pws);
            resultsGOF.gofWholeFit = Obj.OptimizedVPop.gof;
            
            % Run the cross-validation:
            linearCalibrationObjectCrossValidation = Obj.runCrossValidation();
            resultsRMSE.rmseTraining = linearCalibrationObjectCrossValidation.CrossValidationResults.trainingRMSE;
            resultsRMSE.rmseTest = linearCalibrationObjectCrossValidation.CrossValidationResults.testRMSE;
            
            % for debugging:
%             if results.rmse.test < results.rmse.training
%                 nan;
%             end
            
            % Predict against the combo, using the whole-fit calibration
            if ~isempty(p.Results.vPopToPredictAgainst)
                outputDirPrediction = outputDir;
                resultsRMSE.rmsePrediction = Obj.predictAgainstSpecifiedVPop(p.Results.vPopToPredictAgainst,'saveFigDir',outputDirPrediction,'closeFigsAutomatically',true,'figFileNamePrefix','predictionVPop');
            end

            % Save the objects:
            save([outputDir 'linearCalibrationObject.mat'],'Obj');
            save([outputDir 'linearCalibrationObjectCrossValidation.mat'],'linearCalibrationObjectCrossValidation');

            % Generate plots of the fits for this fraction:
            Obj.plotInputVPopFits('saveFigDir',outputDir,'closeFigsAutomatically',true);
            Obj.plotOptimizedVPopFits('saveFigDir',outputDir,'closeFigsAutomatically',true);
            Obj.plotScatterOfOptimVsInputPrevalenceWeights('saveFigDir',outputDir,'closeFigsAutomatically',true);
            
        end
    
    end
    
    function rmse = predictAgainstSpecifiedVPop(Obj,myVPop,varargin)
    % Calculates the root-mean-squared error (RMSE) for the model's prediction of
    % the data specified in a VPop. This function does not perform any fit;
    % it only assesses its prediction according to the optimal prevalence
    % weights. In addition to calculating the RMSE, this function will also
    % generate and save figures of the prediction.
    %
    % INPUT ARGUMENTS
    % -myVPop (VPop object): Contains the data against which the
    % predictions will be tested.
    % -saveFigDir: (string, optional) Directory where to save figures. Does not
    % include the filename (the filename is generated automatically). But
    % it must include a slash at the end of the string. Default is '',
    % which generates the figure in the working directory. Specified as a
    % name-value pair.
    % -closeFigsAutomatically: (boolean, optional) Specifies whether the figures
    % should be automatically closed after they are generated. Default is
    % 'false'. Specified as a name-value pair.
    % -figFileNamePrefix: (character vector) prefix of saved figure
    % filename. Default is 'optimizedVPop'.
    % 
    
        % Our approach will be to create a new object ('ObjNewVPop') corresponding to
        % 'myVPop' and construct the matrices. Intead of fitting the
        % constructed matrices, we set the prevalence weights equal to the
        % ones in 'Obj'. Then, we calculate the RMSE and produce figures of
        % the "fit" of the prediction.
        
        % Instantiate an object corresponding to the specified VPop:
        ObjNewVPop = LinearCalibration(myVPop,'optimOptions',Obj.OptimOptions);
        
        % Construct the matrices:
        ObjNewVPop = ObjNewVPop.constructLinearProblemMatrices();
        
        % Bypass the optimization:
        ObjNewVPop.OptimalPrevalenceWeightsNormalized = Obj.OptimalPrevalenceWeightsNormalized;
        
        % Calculate the RMSE:
        [~,rmse,~] = ObjNewVPop.calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs();
        
        % Below, we use function that were initially developped to be run
        % post-optimization, but we are using them here to assess the
        % prediction (since we bypassed the optimization above):
        ObjNewVPop = ObjNewVPop.constructNewVPopPostOptimization();
        ObjNewVPop.plotOptimizedVPopFits(varargin{:});
    end
    
    function [residuals,rmse,weights] = evaluateDataGroupResponse(Obj,dataGroup)
    % Calculates the residuals and RMSE for the specified dataGroup
    %
    % INPUT ARGUMENTS
    % -dataGroup: (scalar)
    %
        
        % Isolate the specified data group from the linear problem
        % matrices, and then calculate the statistics:
        ObjParticular = Obj.selectFromLinearProblemMatrices('includedDataGroups',dataGroup);
        [residuals,rmse,weights] = ...
            ObjParticular.calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs();
    end
    
    
    function Obj = calculateInputVPopStatistics(Obj,varargin)
    % Computes the residuals and root-mean-squared error for the input VPop
    % prevalence weights
    
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'vPopToPredictAgainst',[]);
        addOptional(p,'saveFigDir','inputVPopPredictions/');
        addOptional(p,'closeFigsAutomatically',false);
        parse(p,varargin{:});
    
        % Instantiate an object corresponding to the Input VPop:
        ObjInputVPop = LinearCalibration(Obj.InputVPop,'optimOptions',Obj.OptimOptions);
        
        ObjInputVPop.OptimalPrevalenceWeightsNormalized = Obj.InputVPop.pws;
        
        % Construct the matrices:
        if ~isempty(Obj.LinearProblemMatrices)
            ObjInputVPop.LinearProblemMatrices = Obj.LinearProblemMatrices;
        else
            ObjInputVPop = ObjInputVPop.constructLinearProblemMatrices();
        end
        
        % calculate residuals and RMSE:
        [Obj.InputVPopStatistics.residuals,Obj.InputVPopStatistics.rootMeanSquaredError,~] = ...
            ObjInputVPop.calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs();
        
        % If indicated, perform prediction against the specified VPop:
        if ~isempty(p.Results.vPopToPredictAgainst)
            Obj.InputVPopStatistics.predictionRootMeanSquaredError = ...
                ObjInputVPop.predictAgainstSpecifiedVPop(p.Results.vPopToPredictAgainst,'saveFigDir',p.Results.saveFigDir,'closeFigsAutomatically',p.Results.closeFigsAutomatically,'figFileNamePrefix','predictionVPop');
        end
    end
    
    function plotInputVPopFits(Obj,varargin)
    % Figures of input VPop fits
    %
    % INPUT ARGUMENTS
    % -saveFigDir: (string, optional) Directory where to save figures. Does not
    % include the filename (the filename is generated automatically). But
    % it must include a slash at the end of the string. Default is '',
    % which generates the figure in the working directory. Specified as a
    % name-value pair.
    % -closeFigsAutomatically: (boolean, optional) Specifies whether the figures
    % should be automatically closed after they are generated. Default is
    % 'false'. Specified as a name-value pair.
    % -figFileNamePrefix: (character vector) prefix of saved figure
    % filename. Default is 'inputVPop'.
    %
    
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'saveFigDir','');
        addOptional(p,'closeFigsAutomatically',false);
        addOptional(p,'figFileNamePrefix','inputVPop');
        parse(p,varargin{:});

        figFilepathPrefix = [p.Results.saveFigDir p.Results.figFileNamePrefix];
        if ~isempty(p.Results.saveFigDir) && ~exist(p.Results.saveFigDir,'dir')
            mkdir(p.Results.saveFigDir);
        end
        createVPopPlots(Obj.InputVPop,p.Results.figFileNamePrefix,figFilepathPrefix,p.Results.closeFigsAutomatically)
    end
    
    function plotOptimizedVPopFits(Obj,varargin)
    % Figures of new (optimized) VPop fits
    %
    % INPUT ARGUMENTS
    % -saveFigDir: (string, optional) Directory where to save figures. Does not
    % include the filename (the filename is generated automatically). But
    % it must include a slash at the end of the string. Default is '',
    % which generates the figure in the working directory. Specified as a
    % name-value pair.
    % -closeFigsAutomatically: (boolean, optional) Specifies whether the figures
    % should be automatically closed after they are generated. Default is
    % 'false'. Specified as a name-value pair.
    % -figFileNamePrefix: (character vector) prefix of saved figure
    % filename. Default is 'optimizedVPop'.
    %
    
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'saveFigDir','');
        addOptional(p,'closeFigsAutomatically',false);
        addOptional(p,'figFileNamePrefix','optimizedVPop');
        parse(p,varargin{:});
    
        figFilepathPrefix = [p.Results.saveFigDir p.Results.figFileNamePrefix];
        if ~isempty(p.Results.saveFigDir) && ~exist(p.Results.saveFigDir,'dir')
            mkdir(p.Results.saveFigDir);
        end
        createVPopPlots(Obj.OptimizedVPop,p.Results.figFileNamePrefix,figFilepathPrefix,p.Results.closeFigsAutomatically)
    end
    
    function plotScatterOfOptimVsInputPrevalenceWeights(Obj,varargin)
	% Scatter plot of new vs old prevalence weights
    %
    % INPUT ARGUMENTS
    % -saveFigDir: (string, optional) Directory where to save figures. Does not
    % include the filename (the filename is generated automatically). But
    % it must include a slash at the end of the string. Default is '',
    % which generates the figure in the working directory. Specified as a
    % name-value pair.
    % -closeFigsAutomatically: (boolean, optional) Specifies whether the figures
    % should be automatically closed after they are generated. Default is
    % 'false'. Specified as a name-value pair.
    % -figFileNamePrefix: (character vector) prefix of saved figure
    % filename. Default is 'prevalenceWeightsNewVsOldScatterPlot'.
    %
        
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'saveFigDir','');
        addOptional(p,'closeFigsAutomatically',false);
        addOptional(p,'figFileNamePrefix','prevalenceWeightsNewVsOldScatterPlot');
        parse(p,varargin{:});
    
        if ~isempty(p.Results.saveFigDir) && ~exist(p.Results.saveFigDir,'dir')
            mkdir(p.Results.saveFigDir);
        end
        figHandle = figure('Name','prevalence weight optimal vs input');
        scatter(Obj.InputVPop.pws,Obj.OptimizedVPop.pws,'LineWidth',1.5);
        xlabel('prevalence weights in input VPop');
        ylabel('optimal prevalence weights');
        ax = gca();
        ax.FontSize = 20;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        linearCalibrationSaveFig(figHandle,[p.Results.saveFigDir p.Results.figFileNamePrefix]);
        if p.Results.closeFigsAutomatically
            close(figHandle);
        end
    end
    
    function plotPrevalenceWeightConfidenceIntervals(Obj,varargin)
    % Plots the confidence intervals of the prevalence weights.
    % 
    % INPUT ARGUMENTS
    % -saveFigDir: (string, optional) Directory where to save figures. Does not
    % include the filename (the filename is generated automatically). But
    % it must include a slash at the end of the string. Default is '',
    % which generates the figure in the working directory. Specified as a
    % name-value pair.
    % -closeFigsAutomatically: (boolean, optional) Specifies whether the figures
    % should be automatically closed after they are generated. Default is
    % 'false'. Specified as a name-value pair.
    
        % Parse the optional arguments:
        p = inputParser;
        addOptional(p,'saveFigDir','');
        addOptional(p,'closeFigsAutomatically',false);
        parse(p,varargin{:});
    
        % Make plot of confidence intervals:
        if ~isempty(p.Results.saveFigDir) && ~exist(p.Results.saveFigDir,'dir')
            mkdir(p.Results.saveFigDir);
        end
        vpIndices = 1:Obj.NumVPs;
        figHandle = figure;
        lineObjects = plot([vpIndices;vpIndices],Obj.OptimalPrevalenceWeightsNormalizedConfidenceIntervals);
        lineColors = [lineObjects(:).Color];
        lineColors = (reshape(lineColors,[3 Obj.NumVPs]))';
        hold on;
        scatter(vpIndices,Obj.OptimalPrevalenceWeightsNormalized,[],lineColors);
        ax = gca();
        ax.FontSize = 20;
%         ax.YScale = 'log';
        ax.YGrid = 'on';
        xlabel('VP index');
        ylabel('optimal prevalence weight and 95% CI');
        xlimPadding = 40;
        xlim([-xlimPadding Obj.NumVPs+xlimPadding]);
        linearCalibrationSaveFig(figHandle,[p.Results.saveFigDir 'prevalenceWeightConfidenceIntervals']);
        if p.Results.closeFigsAutomatically
            close(figHandle);
        end
        
        figHandle = figure;
        ciWidth = Obj.OptimalPrevalenceWeightsNormalizedConfidenceIntervals(2,:) - Obj.OptimalPrevalenceWeightsNormalizedConfidenceIntervals(1,:);
        scatter(Obj.OptimalPrevalenceWeightsNormalized',ciWidth);
        ax = gca();
        ax.FontSize = 20;
        ax.XScale = 'log';
%         ax.YScale = 'log';
        ax.YGrid = 'on';
        xlabel('optimal prevalence weight');
        ylabel('width of 95% CI');
        linearCalibrationSaveFig(figHandle,[p.Results.saveFigDir 'prevalenceWeightConfidenceIntervalWidthsVsEstimate']);
        if p.Results.closeFigsAutomatically
            close(figHandle);
        end
    end

end


%% STATIC METHODS

methods (Static)
    function effN = calculateEffN(prevalenceWeights)
        if abs(sum(prevalenceWeights)-1) > 1e-12
            error('Prevalence weights must sum to 1.');
        end
        effN = 2^(-nansum(prevalenceWeights.*log2(prevalenceWeights)));
    end
end

%% HIDDEN METHODS

methods (Hidden = true)
    
    function Obj = constructLinearProblemMatrices(Obj)
    % Constructs the matrices for the linearized fitting problem.

        % The prevalence weights are required to renormalize VP prevalence weights
        % in cases where some VP simulation values are NaN since they dropped off
        % therapy. With VPs missing, the prevalence weights need to be renormalized
        % before calculating the weighted mean and standard deviation of a
        % simulation variable across the virtual population.
        ignoreNaN = false;
        if strcmpi(Obj.OptimOptions.priorPrevalenceWeightAssumption,"uniform")
            pws = (1/Obj.NumVPs)*ones(1,Obj.NumVPs);
        elseif strcmpi(Obj.OptimOptions.priorPrevalenceWeightAssumption,"specified")
            pws = Obj.InputVPop.pws;
        elseif strcmpi(Obj.OptimOptions.priorPrevalenceWeightAssumption,"ignoreDropout")
            pws = NaN;
            ignoreNaN = true;
        else
            error('The specified optimOptions.priorPrevalenceWeightAssumption is not supported.');
        end

        % Convert the VPop simData into a table, for easier indexing later
        % on:
        simDataRowInfoTbl = cell2table(Obj.InputVPop.simData.rowInfo,'VariableNames',Obj.InputVPop.simData.rowInfoNames);

        % Define a struct array of length 'nDataGroups', where a data group
        % contains all of the datapoints that are dependent (not
        % independent) on each other (e.g., all of the datapoints
        % corresponding to a CDF). Each element of this array is a 
        dataConsolidated = [];

        % Extract data for 'binTable'
        if istable(Obj.InputVPop.binTable)
            for iBinTableRow = 1:height(Obj.InputVPop.binTable)
            % Loop through rows in table

                % Determine the index of the row in the simData table that corresponds
                % to this particular 'binTable' row:
                iMaskSimDataRow = ...
                    simDataRowInfoTbl.time == Obj.InputVPop.binTable.time(iBinTableRow) & ...
                    strcmp(simDataRowInfoTbl.interventionID,Obj.InputVPop.binTable.interventionID{iBinTableRow}) & ...
                    strcmp(simDataRowInfoTbl.elementID,Obj.InputVPop.binTable.elementID{iBinTableRow}) & ...
                    strcmp(simDataRowInfoTbl.elementType,Obj.InputVPop.binTable.elementType{iBinTableRow}) & ...
                    strcmp(simDataRowInfoTbl.expVarID,Obj.InputVPop.binTable.expVarID{iBinTableRow});

                % extract bin edges
                binTableVariableNames = Obj.InputVPop.binTable.Properties.VariableNames;
				binEdges = [-Inf, Obj.InputVPop.binTable{iBinTableRow,'binEdges'}{1}, Inf];

                % Extract the simulation values, and bin them according to the binEdges
                % determined above
                simVals = Obj.InputVPop.simData.Data(iMaskSimDataRow,:);
                simValsBinned = arrayfun(@(x) find(binEdges<x,1,'last'),simVals);

                nObservations = length(binEdges)-1;
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.binTableGroupWeight);
                % Loop through bins and define an observation row for each bin
                for iBin = 1:nObservations
                    % The observed response is taken to be the experimentally observed
                    % probability of this bin
					observationValParticularBin = Obj.InputVPop.binTable{iBinTableRow,'expBins'}{1};
					observationValParticularBin = observationValParticularBin(iBin);
                    % Description of this observation:
                    descriptionParticular = ['binTable; Row ' num2str(iBinTableRow) '; Bin ' num2str(iBin)];
                    % The values for the independent variables will be set to 1 for
                    % VPs in the bin, and 0 for those outside the bin:
                    independentVarValsParticularBin = simValsBinned;
                    independentVarValsParticularBin(independentVarValsParticularBin~=iBin) = 0;
                    independentVarValsParticularBin(independentVarValsParticularBin==iBin) = 1;
                    dataParticular.independentVarVals(iBin,:) = independentVarValsParticularBin;
                    dataParticular.observationVals(iBin) = observationValParticularBin;
                    dataParticular.observationWeights(iBin) = Obj.InputVPop.binTable.weight(iBinTableRow);
                    dataParticular.observationDescriptions{iBin} = descriptionParticular;
                    dataParticular.expWeight(iBin) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.binTable.expN(iBinTableRow),nan,descriptionParticular);
                end
                dataConsolidated = [dataConsolidated dataParticular];
            end
        end

        % Extract data from distTable
        if istable(Obj.InputVPop.distTable)
            nRows = height(Obj.InputVPop.distTable);
            for iDistTableRow = 1:nRows
            % loop through table rows
                % Some of the code below is adapted from the QSP Toolbox function
                % 'plotDistCDFVPop.m':
                myTable = Obj.InputVPop.distTable;
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
                if isstring(Obj.OptimOptions.cdfProbsToFit) && strcmpi(Obj.OptimOptions.cdfProbsToFit,"all")
                    cdfProbsToFit = CDFexp;
                elseif isnumeric(Obj.OptimOptions.cdfProbsToFit)
                    cdfProbsToFit = Obj.OptimOptions.cdfProbsToFit;
                else
                    error('Unsupported specification for Obj.OptimOptions.cdfProbsToFit');
                end
                nObservations = length(cdfProbsToFit);
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.distTableGroupWeight);
                for iIncludedProbs = 1:nObservations
                % loop through the probabilities to fit
                    targetProb = cdfProbsToFit(iIncludedProbs);
                    % find the index that corresponds to the closest probability to the
                    % target probability:
                    closestIndicesToTargetProb = find(abs(CDFexp-targetProb) == min(abs(CDFexp-targetProb)));
                    % In case there are multiple closest indices, take the middle of
                    % the indices:
                    middleClosestIndexToTargetProb = round(closestIndicesToTargetProb(1) + (closestIndicesToTargetProb(end)-closestIndicesToTargetProb(1))/2);
                    % The response value is the experimentally observed CDF
                    % probability:
                    observationValParticular = CDFexp(middleClosestIndexToTargetProb);
                    descriptionParticular = ['distTable; Row ' num2str(iDistTableRow) '; Included Prob ' num2str(iIncludedProbs)];
                    % The values of the independent variables will be set to 1 for
                    % VPs that contribute to the particular point on the CDF, and
                    % to 0 for VPs that don't contribute to that point:
                    if isnan(Obj.InputVPop.distTable.simCombinedIndices{iDistTableRow}(middleClosestIndexToTargetProb))
                        vpIndicesInBin = [];
                    else
                        vpIndicesInBin = Obj.InputVPop.distTable.predIndices{iDistTableRow}(1:Obj.InputVPop.distTable.simCombinedIndices{iDistTableRow}(middleClosestIndexToTargetProb));
                    end
                    independentVarValParticular = zeros(1,Obj.NumVPs);
                    independentVarValParticular(vpIndicesInBin) = 1;
                    dataParticular.independentVarVals(iIncludedProbs,:) = independentVarValParticular;
                    dataParticular.observationVals(iIncludedProbs) = observationValParticular;
                    dataParticular.observationWeights(iIncludedProbs) = Obj.InputVPop.distTable.weight(iDistTableRow);
                    dataParticular.observationDescriptions{iIncludedProbs} = descriptionParticular;
                    dataParticular.expWeight(iIncludedProbs) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.distTable.expN(iDistTableRow),nan,descriptionParticular);
                end
                dataConsolidated = [dataConsolidated dataParticular];
            end
        end

        % Extract data related to the 'brTable'
        % RECIST Class definitions
        % 0 = CR, 1 = PR, 2 = SD, 3 = PD
        recistClassStr = {'CR','PR','SD','PD'};
        resistClassNum = [0 1 2 3];
        if istable(Obj.InputVPop.brTableRECIST)
            for iBRTblRow = 1:height(Obj.InputVPop.brTableRECIST)
            % Loop through table rows
                % Extract the simulation data:
                iMaskSimDataRow = Obj.InputVPop.simData.brRows == iBRTblRow;
                brDataParticular = Obj.InputVPop.simData.brData(iMaskSimDataRow,:);
                nObservations = length(resistClassNum);
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.brTableRECISTGroupWeight);
                for iRECISTClass = 1:nObservations
                % loop through RECIST classes
                    % The response value is the experimentally observed probability for
                    % this RECIST class:
                    observationValParticular = Obj.InputVPop.brTableRECIST.(['exp' recistClassStr{iRECISTClass}])(iBRTblRow);
                    descriptionParticular = ['brTableRECIST; Row ' num2str(iBRTblRow) '; RECIST Class ' num2str(iRECISTClass)];
                    % The values of the independent variables will be set to 1 for
                    % VPs that are in the RECIST class, and 0 for those that aren't
                    vpsInClass = brDataParticular == resistClassNum(iRECISTClass);
                    dataParticular.independentVarVals(iRECISTClass,:) = vpsInClass;
                    dataParticular.observationVals(iRECISTClass) = observationValParticular;
                    dataParticular.observationWeights(iRECISTClass) = Obj.InputVPop.brTableRECIST.weight(iBRTblRow);
                    dataParticular.observationDescriptions{iRECISTClass} = descriptionParticular;
                    dataParticular.expWeight(iRECISTClass) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.brTableRECIST.expN(iBRTblRow),nan,descriptionParticular);
                end
                dataConsolidated = [dataConsolidated dataParticular];
            end
        end
        
        if istable(Obj.InputVPop.rTableRECIST)
            for iRTblRow = 1:height(Obj.InputVPop.rTableRECIST)
                iMaskSimDataRow = Obj.InputVPop.simData.rRows == iRTblRow;
                rDataParticular = Obj.InputVPop.simData.rData(iMaskSimDataRow,:);
                nObservations = length(resistClassNum);
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.rTableRECISTGroupWeight);
                for iRECISTClass = 1:nObservations
                % loop through RECIST classes
                    % The response value is the experimentally observed probability for
                    % this RECIST class:
                    observationValParticular = Obj.InputVPop.rTableRECIST.(['exp' recistClassStr{iRECISTClass}])(iRTblRow);
                    descriptionParticular = ['rTableRECIST; Row ' num2str(iRTblRow) '; RECIST Class ' num2str(iRECISTClass)];
                    % The values of the independent variables will be set to 1 for
                    % VPs that are in the RECIST class, and 0 for those that aren't
                    vpsInClass = rDataParticular == resistClassNum(iRECISTClass);
                    dataParticular.independentVarVals(iRECISTClass,:) = vpsInClass;
                    dataParticular.observationVals(iRECISTClass) = observationValParticular;
                    dataParticular.observationWeights(iRECISTClass) = Obj.InputVPop.rTableRECIST.weight(iRTblRow);
                    dataParticular.observationDescriptions{iRECISTClass} = descriptionParticular;
                    dataParticular.expWeight(iRECISTClass) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.rTableRECIST.expN(iRTblRow),nan,descriptionParticular);
                end
                dataConsolidated = [dataConsolidated dataParticular];
            end
        end

        % Extract data for mean and standard deviation
        if istable(Obj.InputVPop.mnSDTable)
            for iSumStatRow = 1:height(Obj.InputVPop.mnSDTable)
            % Loop through table rows
                nObservations = 2;
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.mnSDTableGroupWeight);

                % Extract the appropriate simulation data:
                iMaskSimDataRow = ...
                    simDataRowInfoTbl.time == Obj.InputVPop.mnSDTable.time(iSumStatRow) & ...
                    strcmp(simDataRowInfoTbl.interventionID,Obj.InputVPop.mnSDTable.interventionID{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementID,Obj.InputVPop.mnSDTable.elementID{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementType,Obj.InputVPop.mnSDTable.elementType{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.expVarID,Obj.InputVPop.mnSDTable.expVarID{iSumStatRow});
                simVals = Obj.InputVPop.simData.Data(iMaskSimDataRow,:);

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
                dataParticular.observationVals(1) = Obj.InputVPop.mnSDTable.expMean(iSumStatRow);
                dataParticular.observationDescriptions{1} = ['mnSDTable; Row ' num2str(iSumStatRow) '; mean'];
                % The values for the independent variables are simply the
                % simulation values
                dataParticular.independentVarVals(1,:) = simVals;
                % renormalize prevalence weights wtih NaN VPs purged. We will do
                % this by scaling the values for the independent variables:
                dataParticular.independentVarVals(1,:) = (1/sumPWsForNonNanVPs)*dataParticular.independentVarVals(1,:); 
                % NaN values will cause issues in the optimization algorithm. 
                % Since we renormalized the prevalence weight, change NaN values to
                % zero (which shouldn't affect the calculation):
                dataParticular.independentVarVals(1,nanValsIndexMask) = 0;
                dataParticular.observationWeights(1) = Obj.InputVPop.mnSDTable.weightMean(iSumStatRow);
                dataParticular.expWeight(1) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.mnSDTable.expN(iSumStatRow),Obj.InputVPop.mnSDTable.expSD(iSumStatRow),dataParticular.observationDescriptions{1});

                % Fit SD (actually, variance)
                % The response value will be the variance:
                dataParticular.observationVals(2) = Obj.InputVPop.mnSDTable.expSD(iSumStatRow)^2;
                dataParticular.observationDescriptions{2} = ['mnSDTable; Row ' num2str(iSumStatRow) '; variance'];
                % The values for the independent variables are the individual
                % variances (i.e., 'residualsSquared'). We cannot yet know the
                % predicted mean, since that would require prior knowledge of the
                % optimal prevalence weights. But, since we are fitting the mean,

                % let's assume that the predicted mean will be (approximately)
                % equal to the experimental mean.
                expMean = Obj.InputVPop.mnSDTable.expMean(iSumStatRow);
                residualsSquared = (simVals-expMean).^2;
                dataParticular.independentVarVals(2,:) = residualsSquared;
                % renormalize prevalence weights wtih NaN VPs purged. We will do
                % this by scaling the values for the independent variables:
                dataParticular.independentVarVals(2,:) = (1/sumPWsForNonNanVPs)*dataParticular.independentVarVals(2,:); 
                % NaN values will cause issues in the optimization algorithm. 
                % Since we renormalized the prevalence weight, change NaN values to
                % zero (which shouldn't affect the calculation):
                dataParticular.independentVarVals(2,nanValsIndexMask) = 0;
                dataParticular.observationWeights(2) = Obj.InputVPop.mnSDTable.weightSD(iSumStatRow);
                dataParticular.expWeight(2) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.mnSDTable.expN(iSumStatRow),Obj.InputVPop.mnSDTable.expSD(iSumStatRow),dataParticular.observationDescriptions{2});

                dataConsolidated = [dataConsolidated dataParticular];
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract data for 2D correlation
        if istable(Obj.InputVPop.corTable)
            for iSumStatRow = 1:height(Obj.InputVPop.corTable)
            % Loop through table rows
                nObservations = 1;
                dataParticular = initDataParticular(Obj.NumVPs,nObservations,Obj.OptimOptions.corTableGroupWeight);


                % Extract the appropriate simulation data:
                jMaskSimData1Row = ...
                    simDataRowInfoTbl.time == Obj.InputVPop.corTable.time1(iSumStatRow) & ...
                    strcmp(simDataRowInfoTbl.interventionID,Obj.InputVPop.corTable.interventionID1{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementID,Obj.InputVPop.corTable.elementID1{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementType,Obj.InputVPop.corTable.elementType1{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.expVarID,Obj.InputVPop.corTable.expVarID1{iSumStatRow});
                simVals1 = Obj.InputVPop.simData.Data(jMaskSimData1Row,:);

                jMaskSimData2Row = ...
                    simDataRowInfoTbl.time == Obj.InputVPop.corTable.time2(iSumStatRow) & ...
                    strcmp(simDataRowInfoTbl.interventionID,Obj.InputVPop.corTable.interventionID2{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementID,Obj.InputVPop.corTable.elementID2{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.elementType,Obj.InputVPop.corTable.elementType2{iSumStatRow}) & ...
                    strcmp(simDataRowInfoTbl.expVarID,Obj.InputVPop.corTable.expVarID2{iSumStatRow});
                simVals2 = Obj.InputVPop.simData.Data(jMaskSimData2Row,:);


                % Now we need to find the experimental mean for the 2d sample,
                % which are not available in the corTable
                iMaskSample1Row = (Obj.InputVPop.mnSDTable.time== Obj.InputVPop.corTable.time1(iSumStatRow)) & ...
                    strcmp(Obj.InputVPop.mnSDTable.interventionID,Obj.InputVPop.corTable.interventionID1{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.elementID,Obj.InputVPop.corTable.elementID1{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.elementType,Obj.InputVPop.corTable.elementType1{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.expVarID,Obj.InputVPop.corTable.expVarID1{iSumStatRow});
                iMaskSample2Row = (Obj.InputVPop.mnSDTable.time== Obj.InputVPop.corTable.time2(iSumStatRow)) & ...
                    strcmp(Obj.InputVPop.mnSDTable.interventionID,Obj.InputVPop.corTable.interventionID2{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.elementID,Obj.InputVPop.corTable.elementID2{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.elementType,Obj.InputVPop.corTable.elementType2{iSumStatRow}) & ...
                    strcmp(Obj.InputVPop.mnSDTable.expVarID,Obj.InputVPop.corTable.expVarID2{iSumStatRow});
                % now we can extract exp mean for sample 1 and sample 2
                sampleexpmean = [Obj.InputVPop.mnSDTable.expMean(iMaskSample1Row);Obj.InputVPop.mnSDTable.expMean(iMaskSample2Row)];
                sampleexpsd = [Obj.InputVPop.mnSDTable.expSD(iMaskSample1Row);Obj.InputVPop.mnSDTable.expSD(iMaskSample2Row)];
                % calculate the approximated correlation that will go into "A";
                % note that here we replaced sample mean/sd with exp mean/sd in
                % order to translate into a linear problem
                simValscorr = ((simVals1-sampleexpmean(1)).*(simVals2-sampleexpmean(2)))/(sampleexpsd(1)*(sampleexpsd(2)));



                % Some simVals will be NaN for VPs that dropped off of therapy. Thus,
                % we need to scale the prevalence weights so that the sum of the
                % weights for the VPs still on therapy adds to 1.
                nanValsIndexMask = isnan(simVals1) | isnan(simVals2);
                if ignoreNaN
                    if any(nanValsIndexMask)
                        continue;
                    end
                    % If we reached here, there are no NaN VPs and thus sumPWsForNonNanVPs = 1
                    sumPWsForNonNanVPs = 1;
                else
                    sumPWsForNonNanVPs = sum(pws(~nanValsIndexMask));
                end

                % Fit the correlations
                % The response value is the experimental correlation:
                dataParticular.observationVals(1) = Obj.InputVPop.corTable.expCor(iSumStatRow);
                dataParticular.observationDescriptions{1} = ['corTable; Row ' num2str(iSumStatRow) '; correlation'];
                % The values for the independent variables are simply the
                % simulation values
                dataParticular.independentVarVals(1,:) = simValscorr;
                % renormalize prevalence weights wtih NaN VPs purged. We will do
                % this by scaling the values for the independent variables:
                dataParticular.independentVarVals(1,:) = (1/sumPWsForNonNanVPs)*dataParticular.independentVarVals(1,:); 
                % NaN values will cause issues in the optimization algorithm. 
                % Since we renormalized the prevalence weight, change NaN values to
                % zero (which shouldn't affect the calculation):
                dataParticular.independentVarVals(1,nanValsIndexMask) = 0;
                dataParticular.observationWeights(1) = Obj.InputVPop.corTable.weight(iSumStatRow);
                dataParticular.expWeight(1) = Obj.OptimOptions.expWeightFuncHandle(Obj.InputVPop.corTable.expN(iSumStatRow),nan,dataParticular.observationDescriptions{1});
                dataConsolidated = [dataConsolidated dataParticular];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
        
        
        

        % Loop through the data groups and perform the following
        % operations:
        % 1. Ensure that there are no missing values.
        % 2. Delete observations which are zero if the objective is
        % calculated relative to the observation value (since we don't want
        % to divide by zero)
        % 3. Calculate the number of observation in the data group.
        Obj.IgnoredObservationsDescriptions = {};
        numObservationsInEachDataGroup = nan(1,length(dataConsolidated));
        for iDataCons = 1:length(dataConsolidated)
            % As a sanity check, ensure that none of the fields in
            % 'dataConsolidated have any missing values:
            fieldNames = fields(dataConsolidated(iDataCons));
            for iField = 1:length(fieldNames)
                if ( iscell(dataConsolidated(iDataCons).(fieldNames{iField})) && any(cellfun(@isempty,dataConsolidated(iDataCons).(fieldNames{iField})),'all') ) || ...
                   ( isnumeric(dataConsolidated(iDataCons).(fieldNames{iField})) && any(isnan(dataConsolidated(iDataCons).(fieldNames{iField})),'all') )
                    error('Missing data values were produced when constructing the linear problem matrices.');
                end
            end

            % delete observations = 0 if fit is relative
            if strcmpi(Obj.OptimOptions.responseValTransformation,"relative")
                iMaskZeroObs = dataConsolidated(iDataCons).observationVals == 0;
                Obj.IgnoredObservationsDescriptions = ...
                    [Obj.IgnoredObservationsDescriptions; dataConsolidated(iDataCons).observationDescriptions(iMaskZeroObs)];
                dataConsolidated(iDataCons).independentVarVals(iMaskZeroObs,:) = [];
                dataConsolidated(iDataCons).observationVals(iMaskZeroObs) = [];
                dataConsolidated(iDataCons).observationWeights(iMaskZeroObs) = [];
                dataConsolidated(iDataCons).observationDescriptions(iMaskZeroObs) = [];
                dataConsolidated(iDataCons).expWeight(iMaskZeroObs) = [];
            end

            numObservationsInEachDataGroup(iDataCons) = length(dataConsolidated(iDataCons).observationVals);
        end

        % In 'dataConsolidated', the data is separated according to data
        % groups. We still need to consolidate this data into large
        % matrices in order to fit the optimization problem:
        
        % The cumulative number of observations across the data groups:
        cumulativeNumObservations = cumsum(numObservationsInEachDataGroup);
        cumulativeNumObservations = [0 cumulativeNumObservations];
        
        % Initialize 'A' in Ax=b
        % Will eventually be matrix of size 'nObservations' x 'Obj.NumVPs'
        Obj.LinearProblemMatrices.independentVarVals = nan(cumulativeNumObservations(end),Obj.NumVPs);
        
        % Initialize 'b' in Ax=b
        % Will eventually be vector of size 'nObservations' x 1
        Obj.LinearProblemMatrices.observationVals = nan(cumulativeNumObservations(end),1);
        
        % Each observation has associated a weight
        % Will eventually be vector of size 'nObservations' x 1
        Obj.LinearProblemMatrices.observationWeights = nan(cumulativeNumObservations(end),1);
        Obj.LinearProblemMatrices.dataGroupWeights = nan(cumulativeNumObservations(end),1);
        
        % Keep track of the data group corresponding to each row:
        Obj.LinearProblemMatrices.dataGroup = nan(cumulativeNumObservations(end),1);
        
        % String descriptions of the different observations -- useful for
        % debugging.
        % Will eventually be vector of size 'nObservations' x 1
        Obj.LinearProblemMatrices.observationDescriptions = cell(cumulativeNumObservations(end),1);
        
        Obj.LinearProblemMatrices.numObservationsInDataGroup = nan(cumulativeNumObservations(end),1);
        
        Obj.LinearProblemMatrices.expWeight = nan(cumulativeNumObservations(end),1);
        
        % Now, construct the large matrices, consolidating all of the data groups:
        for iDataCons = 1:length(dataConsolidated)
            startRowIndex = cumulativeNumObservations(iDataCons)+1;
            endRowIndex = cumulativeNumObservations(iDataCons+1);
            Obj.LinearProblemMatrices.independentVarVals(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).independentVarVals;
            Obj.LinearProblemMatrices.observationVals(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).observationVals;
            Obj.LinearProblemMatrices.observationWeights(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).observationWeights;
            Obj.LinearProblemMatrices.dataGroupWeights(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).dataGroupWeight;
            Obj.LinearProblemMatrices.dataGroup(startRowIndex:endRowIndex,:) = iDataCons;
            Obj.LinearProblemMatrices.observationDescriptions(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).observationDescriptions;
            Obj.LinearProblemMatrices.numObservationsInDataGroup(startRowIndex:endRowIndex,:) = numObservationsInEachDataGroup(iDataCons);
            Obj.LinearProblemMatrices.expWeight(startRowIndex:endRowIndex,:) = dataConsolidated(iDataCons).expWeight;
        end

        % Calculate and incorporate the observation weights:
        Obj.LinearProblemMatrices = incorporateLinearProblemMatricesWeights(Obj.LinearProblemMatrices,Obj.OptimOptions);
        
        % Nested function to initialize a 'dataParticular' struct:
        function dataParticular = initDataParticular(nVPs,nObservations,dataGroupWeight)
            dataParticular.independentVarVals = nan(nObservations,nVPs); % 'A' in Ax=b
            dataParticular.observationVals = nan(nObservations,1); % 'b' in Ax=b
            dataParticular.dataGroupWeight = dataGroupWeight; % weight for the data group
            dataParticular.observationWeights = nan(nObservations,1); % weight for the observation (in addition to the data group)
            dataParticular.observationDescriptions = cell(nObservations,1); % description of the observation
            dataParticular.expWeight = nan(nObservations,1); % experimental sample size
        end
    end
    
    function Obj = selectFromLinearProblemMatrices(Obj,varargin)
    % Selects only particular data groups and VPs from 'Obj.LinearProblemMatrices',
    % and stores the results in 'Obj.LinearProblemMatricesParticular'.
    % Optional input arguments are 'includedDataGroups' and
    % 'includedVPIndices'. Results are stored in
    % Obj.LinearProblemMatricesParticular.
    %
    
        % Parse optional arguments
        in = inputParser();
        in.addParameter('includedDataGroups',[],@isnumeric);
        in.addParameter('includedVPIndices',[],@isnumeric);
        in.parse(varargin{:});
    
        fieldNames = fields(Obj.LinearProblemMatrices);
        nFields = length(fieldNames);
        
        if isempty(in.Results.includedDataGroups)
            Obj.LinearProblemMatricesParticular = Obj.LinearProblemMatrices;
        else
            observationIndicesToInclude = [];
            % Determine which observations to include:
            for iDataGrp = 1:length(in.Results.includedDataGroups)
                observationIndicesToInclude = [observationIndicesToInclude; find(Obj.LinearProblemMatrices.dataGroup == in.Results.includedDataGroups(iDataGrp))];
            end
            % Extract those observations:
            for iField = 1:nFields
                Obj.LinearProblemMatricesParticular.(fieldNames{iField}) = Obj.LinearProblemMatrices.(fieldNames{iField})(observationIndicesToInclude,:); 
            end
        end
        
        if ~isempty(in.Results.includedVPIndices)
            for iField = 1:nFields
                if size(Obj.LinearProblemMatricesParticular.(fieldNames{iField}),2) > 1
                    Obj.LinearProblemMatricesParticular.(fieldNames{iField}) = Obj.LinearProblemMatricesParticular.(fieldNames{iField})(:,in.Results.includedVPIndices);
                end
            end
        end        
    end
    
    function Obj = runOptimization(Obj)
    % Run the optimization. 
    %
        
        % Start timer
        tictoc = tic();
        
        if isempty(Obj.LinearProblemMatricesParticular)
            Obj.LinearProblemMatricesParticular = Obj.LinearProblemMatrices;
        end
        
        nVPs = size(Obj.LinearProblemMatricesParticular.independentVarValsWeighted,2);
        
        % Run the optimization
        if strcmpi(Obj.OptimOptions.optimizationAlgorithm,"lsqnonneg")
            % Set optimization options:
            options = optimset('Display','none');
            if isstruct(Obj.OptimOptions.optimizationAlgorithmOptions)
                fieldNames = fields(Obj.OptimOptions.optimizationAlgorithmOptions);
            else
                fieldNames = {};
            end
            for iField = 1:length(fieldNames)
                options.(fieldNames{iField}) = Obj.OptimOptions.optimizationAlgorithmOptions.(fieldNames{iField});
            end
            % Run the optimization:
            [Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization,Obj.OptimizationResults.normOfresidualsSquared,Obj.OptimizationResults.residuals,Obj.OptimizationResults.exitFlag,Obj.OptimizationResults.output,Obj.OptimizationResults.lagrangeMultipliers] = ...
                lsqnonneg(Obj.LinearProblemMatricesParticular.independentVarValsWeighted,Obj.LinearProblemMatricesParticular.observationValsWeighted,options);
        elseif strcmpi(Obj.OptimOptions.optimizationAlgorithm,"nnls")
            % For poor problems, the function will throw a bunch of
            % warnings, despite still being able to converge. Let's turn
            % off these warnings.
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:rankDeficientMatrix');
            % Run the optimization. Note that this optimization function
            % takes the options as a struct, so we can input them directly.
            [Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization,Obj.OptimizationResults.lagrangeMultipliers,~] = ...
                nnls(Obj.LinearProblemMatricesParticular.independentVarValsWeighted,Obj.LinearProblemMatricesParticular.observationValsWeighted,Obj.OptimOptions.optimizationAlgorithmOptions);
            % turn the warnings back on:
            warning('on','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:rankDeficientMatrix');
            Obj.OptimizationResults.exitFlag = 1;
        elseif strcmpi(Obj.OptimOptions.optimizationAlgorithm,"lsqlin")
            % Specifying an 'x0' would result in a warning that x0 gets ignored in
            % 'lsqlin'
            % x0 = (1/nVPs)*ones(nVPs,1);
            x0 = [];
            lb = zeros(nVPs,1);
            ub = [];
            if ~isempty(Obj.OptimOptions.maxPrevalenceWeight) && isfinite(Obj.OptimOptions.maxPrevalenceWeight)
                A = spdiags(ones(nVPs,1),0,nVPs,nVPs);
                b = Obj.OptimOptions.maxPrevalenceWeight .* ones(1,nVPs);
            else
                A = [];
                b = [];
            end
            % Constraint that sum of prevalence weights should be 1:
            Aeq = ones(1,nVPs);
            beq = 1;
            % Set optimization options:
            lsqlinOptions = optimoptions('lsqlin','Display','none');
            if isstruct(Obj.OptimOptions.optimizationAlgorithmOptions)
                fieldNames = fields(Obj.OptimOptions.optimizationAlgorithmOptions);
            else
                fieldNames = {};
            end
            for iField = 1:length(fieldNames)
                lsqlinOptions.(fieldNames{iField}) = Obj.OptimOptions.optimizationAlgorithmOptions.(fieldNames{iField});
            end
            % Run the optimization:
            [Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization,Obj.OptimizationResults.normOfresidualsSquared,Obj.OptimizationResults.residuals,Obj.OptimizationResults.exitFlag,Obj.OptimizationResults.output,Obj.OptimizationResults.lagrangeMultipliers] = ...
                lsqlin(Obj.LinearProblemMatricesParticular.independentVarValsWeighted,Obj.LinearProblemMatricesParticular.observationValsWeighted,A,b,Aeq,beq,lb,ub,x0,lsqlinOptions);
        else
            error('The specified optimization algorithm is not recognized.');
        end
        
        % if the fit failed to converge, set the solution as NaN:
        if Obj.OptimizationResults.exitFlag < 1
            Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization = nan*Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization;
        end
        
        % normalize the prevalence weights to sum to 1:
        Obj.OptimizationResults.optimalPrevalenceWeightsNormalized = ...
            Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization./sum(Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization);
        
        Obj.OptimizationResults.timeElapsedMinutes = toc(tictoc)/60;
    end
    
    function Obj = performPostOptimizationAnalysis(Obj)
    % Computes some statistics based on the optimal prevalence weight
        
        Obj.OptimizationResults.sumOfOptimalPrevalenceWeightsPriorToNormalization = sum(Obj.OptimizationResults.optimalPrevalenceWeightsPriorToNormalization);
        
        [Obj.OptimizationResults.residualsAfterRenormalizationOfPrevalenceWeights,Obj.OptimizationResults.rmseConsideringOptimalPrevalenceWeightsNormalized,~] = ...
            Obj.calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs();
        
%         % The following calculation follows from the calculations performed
%         % in 'lsqnonneg':
%         optimizationResults.lagrangeMultipliersAfterRenormalizationOfPrevalenceWeights = ...
%             Obj.LinearProblemMatrices.independentVarVals'*optimizationResults.residualsAfterRenormalizationOfPrevalenceWeights;
        
        Obj.OptimizationResults.effectiveN = Obj.calculateEffN(Obj.OptimizationResults.optimalPrevalenceWeightsNormalized);
    end
    
    function Obj = constructNewVPopPostOptimization(Obj)
    % Constructs a new VPop corresponding to the optimal prevalence weights
        
        % Sanity check that prevalence weights sum to one
        if abs(sum(Obj.OptimalPrevalenceWeightsNormalized)-1) > 1e-12
            error('Prevalence weights do not sum to 1 after renormalization.')
        end
        
        % Copy the old VPop and replace the prevalence weights:
        Obj.OptimizedVPop = Obj.InputVPop;
        Obj.OptimizedVPop.pws = (Obj.OptimalPrevalenceWeightsNormalized(:))';

        % recalculate the tables and goodness of fit 
        Obj.OptimizedVPop=addPredTableVals(Obj.OptimizedVPop);
        Obj.OptimizedVPop=evaluateGOF(Obj.OptimizedVPop);
    end
    
    function [residuals,rmse,weights] = calculateResidualsAndRootMeanSquaredErrorForSpecifiedPWs(Obj)
    % Calculate residuals and RMSE corresponding to the specified
    % prevalence weights
        
        % The statistics are calculated according to the matrices in
        % 'Obj.LinearProblemMatricesParticular'. If this property is empty,
        % set it equal to 'Obj.LinearProblemMatrices'
        if isempty(Obj.LinearProblemMatricesParticular)
            Obj.LinearProblemMatricesParticular = Obj.LinearProblemMatrices;
        end
    
        weights = Obj.LinearProblemMatricesParticular.ultimateWeights;
        
        % The actual residuals and RMSE used in the fitting:
        residuals = Obj.LinearProblemMatricesParticular.observationValsWeighted - Obj.LinearProblemMatricesParticular.independentVarValsWeighted*Obj.OptimalPrevalenceWeightsNormalized(:);
        % weighted mean (note that this is a weighted mean since the squared residuals are weighted with weights which sum to 1: (this will be NaN when the weights sum to zero):
        mse = sum(residuals.^2);
        rmse = sqrt(mse);
        
        % for debugging
        % In the fitting, each residual is scaled according to the number
        % of observations in the data group, so that the influence of
        % a data group on the RMSE doesn't depend on the number of
        % observations in the data group. However, if we're analyzing a
        % particular data group's residuals after fitting, it would be more
        % informative to consider the unscaled version, which is calculated
        % below:
%         residualsCorrectedForN = residuals.*sqrt(Obj.LinearProblemMatricesParticular.numObservationsInDataGroup);
%         rmseCorrectedForN = sqrt(mean(residualsCorrectedForN.^2));
        
        % for debugging:
%         uniqueDataGroups = unique(Obj.LinearProblemMatricesParticular.dataGroup);
%         uniqueDataGroups = uniqueDataGroups(:);
%         nDataGroups = length(uniqueDataGroups);
%         rmseCorrectedForNByDataGroup = nan(nDataGroups,1);
%         dataGroupDescriptions = cell(nDataGroups,1);
%         for iGrp = 1:nDataGroups
%             iMaskCrp = Obj.LinearProblemMatricesParticular.dataGroup == iGrp;
%             residualsParticular = residualsCorrectedForN(iMaskCrp);
%             rmseCorrectedForNByDataGroup(iGrp) = sqrt(mean(residualsParticular.^2));
%             observationDescriptionsParticular = Obj.LinearProblemMatricesParticular.observationDescriptions(iMaskCrp);
%             dataGroupDescriptions{iGrp} = observationDescriptionsParticular{1};
%         end
%         rmseByDataGroupTable = table(uniqueDataGroups,dataGroupDescriptions,rmseCorrectedForNByDataGroup,'VariableNames',{'dataGroup','firstObservationDescription','rmseCorrectedForNForDataGroup'});
%         
%         residualsDiagnosticsTable = table(...
%             Obj.LinearProblemMatricesParticular.dataGroup,...
%             Obj.LinearProblemMatricesParticular.observationDescriptions,....
%             residualsCorrectedForN,...
%             Obj.LinearProblemMatricesParticular.numObservationsInDataGroup,...
%             Obj.LinearProblemMatricesParticular.observationVals,...
%             (Obj.LinearProblemMatricesParticular.independentVarVals*prevalenceWeights(:)),...
%             Obj.LinearProblemMatricesParticular.observationValsWeighted.*sqrt(Obj.LinearProblemMatricesParticular.numObservationsInDataGroup),...
%             (Obj.LinearProblemMatricesParticular.independentVarValsWeighted*prevalenceWeights(:)).*sqrt(Obj.LinearProblemMatricesParticular.numObservationsInDataGroup),...
%             'VariableNames',{'dataGroup','observationDescriptions','residualsCorrectedForN','numObservationsInDataGroup','observationVal','predVal','observationValWeighted_x_N','predValWeighted_x_N'});
%         
%         residualsDiagnosticsTable = outerjoin(residualsDiagnosticsTable,rmseByDataGroupTable,'Keys',{'dataGroup'},'RightVariables',{'rmseCorrectedForNForDataGroup'});
%         residualsDiagnosticsTable = sortrows(residualsDiagnosticsTable,'rmseCorrectedForNForDataGroup');
    end

end

end


%% INTERNAL HELPER FUNCTIONS

function createVPopPlots(vPop,figName,figFilepathPrefix,closeFigsAutomatically)
% Internal helper function for creating VPop plots

    plotVPopFunctionHandles = {@plotBinVPop,@plotBRVPop,@plotMnSDVPop,@plotDistCDFVPop};
    for iFunc = 1:length(plotVPopFunctionHandles)
        try
            figHandle = plotVPopFunctionHandles{iFunc}(vPop);
        catch
            warning(['Unable to generate plots using ' func2str(plotVPopFunctionHandles{iFunc}) '. This could be because there is no data of this type to plot.']);
            continue;
        end
        if isempty(figHandle)
            warning(['Unable to generate plots using ' func2str(plotVPopFunctionHandles{iFunc}) '. This could be because there is no data of this type to plot.']);
            continue;
        end
        if ~isempty(figName)
            figHandle.Name = [figName ' - ' func2str(plotVPopFunctionHandles{iFunc})];
        end
        linearCalibrationSaveFig(figHandle,[figFilepathPrefix '_' func2str(plotVPopFunctionHandles{iFunc})]);
        if closeFigsAutomatically
            close(figHandle);
        end
    end
end

function linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension)
% Internal helper function for saving figures
    figure(figHandle);
    set(figHandle, 'Position', [1 1 1680 1050]); % make the figure large before saving
    savefig(figHandle, figFilepathWithoutExtension);
    % export_fig figFilepathWithoutExtension -png -jpg -tif -transparent;
    print(figHandle, [figFilepathWithoutExtension '.tif'],'-dtiff','-r300');
end

function linearProblemMatrices = incorporateLinearProblemMatricesWeights(linearProblemMatrices,optimOptions)
% Internal helper function used during the construction of the linear
% problem matrices for calculating and incorporating the observation
% weights
        
    % Calculate the total weight for each observation:
    % The weights need to be square-rooted, since otherwise they will
    % get squared when the residuals are squared in the least
    % squares routine:
    linearProblemMatrices.ultimateWeights = ...
        sqrt( linearProblemMatrices.expWeight .* linearProblemMatrices.observationWeights .* linearProblemMatrices.dataGroupWeights ./ linearProblemMatrices.numObservationsInDataGroup );

    % normalize so that the squared weights sum to 1
    normalizationDenominator = sqrt(sum(linearProblemMatrices.ultimateWeights.^2));
	% !!! the sqrt can be taken in the following line instead, to be more
	% efficient
    linearProblemMatrices.ultimateWeights = linearProblemMatrices.ultimateWeights/normalizationDenominator;
    
    % Relative transformation of response values, if performed, is performed after 
    % normalization of weights (since we don't want the value of one
    % observation to change the weight of another, and since it is not
    % correct to normalize based on the transformation)
    if strcmpi(optimOptions.responseValTransformation,"relative")
        % I don't think the following needs to be square rooted,
        % since we are assuming that the residuals, not the squared
        % residuals, are proportional to absolute values.
        linearProblemMatrices.ultimateWeights = ...
            abs(linearProblemMatrices.ultimateWeights ./ linearProblemMatrices.observationVals);
    end
    
    % factor in weights
    ultimateWeightsMatrix = spdiags(linearProblemMatrices.ultimateWeights,0,length(linearProblemMatrices.ultimateWeights),length(linearProblemMatrices.ultimateWeights));
    linearProblemMatrices.independentVarValsWeighted = ultimateWeightsMatrix*linearProblemMatrices.independentVarVals;
    linearProblemMatrices.observationValsWeighted = ultimateWeightsMatrix*linearProblemMatrices.observationVals;

end

function generateAndSaveCrossValidationPlots(outputDir,resultsOut,p)
% Internal helper function for generating the cross-validation plots

    %% Define plot visual parameters:
    plotParameters.baggingRMSE.rmseWholeFit.displayName = 'whole-fit RMSE';
    plotParameters.baggingRMSE.rmseWholeFit.lineStyle = '-';
    plotParameters.baggingRMSE.rmseWholeFit.marker = 'o';
    plotParameters.baggingRMSE.rmseWholeFit.color = [0 0 255]/255;
    plotParameters.baggingRMSE.rmseTraining.displayName = 'CV training RMSE';
    plotParameters.baggingRMSE.rmseTraining.lineStyle = ':';
    plotParameters.baggingRMSE.rmseTraining.marker = '+';
    plotParameters.baggingRMSE.rmseTraining.color = [0 147 147]/255;
    plotParameters.baggingRMSE.rmseTest.displayName = 'CV test RMSE';
    plotParameters.baggingRMSE.rmseTest.lineStyle = '--';
    plotParameters.baggingRMSE.rmseTest.marker = 'square';
    plotParameters.baggingRMSE.rmseTest.color = [0 255 255]/255;
    plotParameters.baggingRMSE.rmsePrediction.displayName = 'prediction RMSE';
    plotParameters.baggingRMSE.rmsePrediction.lineStyle = '-.';
    plotParameters.baggingRMSE.rmsePrediction.marker = 'x';
    plotParameters.baggingRMSE.rmsePrediction.color = [round(255/2) 0 255]/255;
    plotParameters.baggingEffN.effNWholeFit.displayName = 'effective n';
    plotParameters.baggingEffN.effNWholeFit.lineStyle = '-';
    plotParameters.baggingEffN.effNWholeFit.marker = 'o';
    plotParameters.baggingEffN.effNWholeFit.color = [255 0 0]/255;
    plotParameters.baggingGOF.gofWholeFit.displayName = 'whole-fit GOF';
    plotParameters.baggingGOF.gofWholeFit.lineStyle = '-';
    plotParameters.baggingGOF.gofWholeFit.marker = 'o';
    plotParameters.baggingGOF.gofWholeFit.color = [0 255 255]/255;
    plotParameters.inputVPopRMSE = plotParameters.baggingRMSE;
    plotParameters.inputVPopEffN = plotParameters.baggingEffN;
    plotParameters.inputVPopGOF = plotParameters.baggingGOF;
    plotParameters.noBaggingRMSE = plotParameters.baggingRMSE;
    plotParameters.noBaggingEffN = plotParameters.baggingEffN;
    plotParameters.noBaggingGOF = plotParameters.baggingGOF;
    
    
    %% Plot the results RMSE and effective N:
    figfilename = 'baggingCrossValidationRMSEAndEffN';
    figHandle = figure;
    xlabel('fraction VPs per bagging sample')
    yyaxis left;
    hold on;
    makeLinePlots(plotParameters.baggingRMSE,resultsOut.baggingRMSE);
    ylabel('root-mean-squared error (RMSE)');
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    yyaxis right;
    hold on;
    makeLinePlots(plotParameters.baggingEffN,resultsOut.baggingEffN);
    ylabel('effective num. virtual patients (n_{eff})');
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    ax = gca();
    ax.XScale = 'log';
    ax.FontSize = 20;
    grid on;
    figFilepathWithoutExtension = [outputDir figfilename 'WithoutLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    legend;
    figFilepathWithoutExtension = [outputDir figfilename 'WithLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);

    %% Make RMSE/EffN bar plot
    
    % First, plot RMSE:
    
    figfilename = 'barRMSEandEffN';
    
    leftAxisDescription = 'baggingRMSE';
    rightAxisDescription = 'baggingEffN';
        
    subgroupTable = constructSubGroupTable(plotParameters,{leftAxisDescription,rightAxisDescription});
    [uniqueSubTypes,uniqueSubTypesIndices,~] = unique(subgroupTable.descriptionsSubType,'stable');
    subgroupLabels = (subgroupTable.displayNames(uniqueSubTypesIndices))';
    subgroupColors = subgroupTable.barColors(uniqueSubTypesIndices);
    
    groupLabels = {};
    barYLeft = [];

    barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'inputVPopRMSE',1)];
    groupLabels = [groupLabels ; {'input VPop'}];
    
    barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'noBaggingRMSE',1)];
    groupLabels = [groupLabels ; {'no bagging'}];
    
    for iFraction = 1:height(resultsOut.baggingRMSE)
        barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'baggingRMSE',iFraction)];
        groupLabels = [groupLabels ; {['bagging with fraction: ' num2str(resultsOut.baggingRMSE.fractionVPsPerBaggingIteration(iFraction))]}];
    end
    
    figHandle = figure;
    yyaxis left;
    barHandleLeft = bar(barYLeft);
    for iBarLeft = 1:length(barHandleLeft)
        barHandleLeft(iBarLeft).FaceColor = subgroupColors{iBarLeft};
    end
    
    ylabel({'root-mean-squared'; 'error (RMSE)'});
    
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    
    % Second, plot EffN:
    
    barYRight = [];
    
    barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'inputVPopEffN',1)];
    
    barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'noBaggingEffN',1)];
    
    for iFraction = 1:height(resultsOut.baggingEffN)
        barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'baggingEffN',iFraction)];
    end
    
    yyaxis right;
    barHandleRight = bar(barYRight);
    for iBarRight = 1:length(barHandleRight)
        barHandleRight(iBarRight).FaceColor = subgroupColors{iBarRight};
    end
    
    ylabel({'effective num. virtual patients';'(n_{eff}) for the whole-fit case'});
    
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    
    ax = gca();
    ax.XTick = 1:length(groupLabels);
    ax.XTickLabel = groupLabels';
    xtickangle(45);
    ax.FontSize = 20;
    ax.YGrid = 'on';
    figFilepathWithoutExtension = [outputDir figfilename 'WithoutLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    legend(subgroupLabels);
    figFilepathWithoutExtension = [outputDir figfilename 'WithLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    

    %% Make GOF/EffN Bar plot
    
    % First, plot RMSE:
    
    figfilename = 'barGOFandEffN';

    leftAxisDescription = 'baggingGOF';
    rightAxisDescription = 'baggingEffN';
    
    subgroupTable = constructSubGroupTable(plotParameters,{leftAxisDescription,rightAxisDescription});
    [uniqueSubTypes,uniqueSubTypesIndices,~] = unique(subgroupTable.descriptionsSubType,'stable');
    subgroupLabels = (subgroupTable.displayNames(uniqueSubTypesIndices))';
    subgroupColors = subgroupTable.barColors(uniqueSubTypesIndices);
    
    groupLabels = {};
    barYLeft = [];
    
    barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'inputVPopGOF',1)];
    groupLabels = [groupLabels ; {'input VPop'}];
    
    barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'noBaggingGOF',1)];
    groupLabels = [groupLabels ; {'no bagging'}];
    
    for iFraction = 1:height(resultsOut.baggingGOF)
        barYLeft = [barYLeft ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'baggingGOF',iFraction)];
        groupLabels = [groupLabels ; {['bagging with fraction: ' num2str(resultsOut.baggingGOF.fractionVPsPerBaggingIteration(iFraction))]}];
    end
    
    figHandle = figure;
    yyaxis left;
    barHandleLeft = bar(barYLeft);
    for iBarLeft = 1:length(barHandleLeft)
        barHandleLeft(iBarLeft).FaceColor = subgroupColors{iBarLeft};
    end
    
    ylabel({'goodness-of-fit (GOF) p-value';'for the whole-fit case'});
    
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    
    % Second, plot EffN:
    
    barYRight = [];
    
    barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'inputVPopEffN',1)];
    
    barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'noBaggingEffN',1)];
    
    for iFraction = 1:height(resultsOut.baggingEffN)
        barYRight = [barYRight ; constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,'baggingEffN',iFraction)];
    end
    
    yyaxis right;
    barHandleRight = bar(barYRight);
    for iBarRight = 1:length(barHandleRight)
        barHandleRight(iBarRight).FaceColor = subgroupColors{iBarRight};
    end
    
    ylabel({'effective num. virtual patients (n_{eff})';'for the whole-fit case'});
    
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    
    ax = gca();
    ax.XTick = 1:length(groupLabels);
    ax.XTickLabel = groupLabels';
    xtickangle(45);
    ax.FontSize = 20;
    ax.YGrid = 'on';
    figFilepathWithoutExtension = [outputDir figfilename 'WithoutLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    legend(subgroupLabels);
    figFilepathWithoutExtension = [outputDir figfilename 'WithLegend'];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    
    
    %%  Plot the GOF
    if p.Results.calculateGOF
        figfilename = 'baggingGOFAndEffN';
        figHandle = figure;
        xlabel('fraction VPs per bagging sample')
        yyaxis left;
        hold on;
        makeLinePlots(plotParameters.baggingGOF,resultsOut.baggingGOF);
        ylabel({'goodness-of-fit (GOF) p-value';'for the whole-fit case'});
        % Place the bottom of the y-axis at the origin:
        yLimMem = ylim();
        ylim([0 yLimMem(2)]);
        yyaxis right;
        hold on;
        makeLinePlots(plotParameters.baggingEffN,resultsOut.baggingEffN);
        ylabel({'effective num. virtual patients';'(n_{eff}) for the whole-fit case'});
        % Place the bottom of the y-axis at the origin:
        yLimMem = ylim();
        ylim([0 yLimMem(2)]);
        ax = gca();
        ax.XScale = 'log';
        ax.FontSize = 20;
        grid on;
        figFilepathWithoutExtension = [outputDir figfilename];
        linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);
    end

    
    %%  Plot the time taken for a single bagging run:
    figfilename = 'timeForSingleBaggingRun';
    figHandle = figure;
    xlabel('fraction VPs per bagging sample')
    hold on;
    plot(resultsOut.baggingRunTime.fractionVPsPerBaggingIteration,resultsOut.baggingRunTime.calibrationTimeForSingleNonCVRunSeconds,'Marker','o','MarkerSize',12,'LineWidth',2,'DisplayName','run time');
    ylabel('run time for bagging [seconds]');
    % Place the bottom of the y-axis at the origin:
    yLimMem = ylim();
    ylim([0 yLimMem(2)]);
    ax = gca();
    ax.XScale = 'log';
    ax.FontSize = 20;
    grid on;
    figFilepathWithoutExtension = [outputDir figfilename];
    linearCalibrationSaveFig(figHandle,figFilepathWithoutExtension);

    
    %% Nested functions
    
    function makeLinePlots(plotParameters,resultsOutTable)
        fieldNamesInternal = fields(plotParameters);
        for iFieldInternal = 1:length(fieldNamesInternal)
            if any(strcmp(resultsOutTable.Properties.VariableNames,fieldNamesInternal{iFieldInternal}))
                plot(resultsOutTable.fractionVPsPerBaggingIteration,resultsOutTable.(fieldNamesInternal{iFieldInternal}),'LineStyle',plotParameters.(fieldNamesInternal{iFieldInternal}).lineStyle,'Marker',plotParameters.(fieldNamesInternal{iFieldInternal}).marker,'Color',plotParameters.(fieldNamesInternal{iFieldInternal}).color,'MarkerSize',12,'LineWidth',1.5,'DisplayName',plotParameters.(fieldNamesInternal{iFieldInternal}).displayName);
            end
        end
    end
    
    function subgroupTable = constructSubGroupTable(plotParameters,axisDescriptions)
        descriptionsMain = {};
        descriptionsSubType = {};
        displayNames = {};
        barColors = {};
        for iDescription = 1:length(axisDescriptions)
            fieldNamesInternal = fields(plotParameters.(axisDescriptions{iDescription}));
            for iFieldInternal = 1:length(fieldNamesInternal)
                descriptionsMain = [descriptionsMain ; axisDescriptions{iDescription}];
                descriptionsSubType =  [descriptionsSubType ; fieldNamesInternal{iFieldInternal}];
                displayNames = [displayNames ; plotParameters.(axisDescriptions{iDescription}).(fieldNamesInternal{iFieldInternal}).displayName];
                barColors = [barColors ; plotParameters.(axisDescriptions{iDescription}).(fieldNamesInternal{iFieldInternal}).color];
            end
        end
        subgroupTable = table(descriptionsMain,descriptionsSubType,displayNames,barColors);
    end
    
    function barYParticularRow = constructBarYRow(resultsOut,subgroupTable,uniqueSubTypes,descriptionMain,iResult)
        iMaskDescriptionMain = strcmp(subgroupTable.descriptionsMain,descriptionMain);
        descriptionSubTypes = subgroupTable.descriptionsSubType(iMaskDescriptionMain);
        barYParticularRow = nan(1,length(uniqueSubTypes));
        
        for iSubType = 1:length(uniqueSubTypes)
            if isfield(resultsOut,descriptionMain) && any(strcmp(resultsOut.(descriptionMain).Properties.VariableNames,uniqueSubTypes{iSubType}))
                barYParticularRow(iSubType) = resultsOut.(descriptionMain).(uniqueSubTypes{iSubType})(iResult);
            end
        end
    end
end
