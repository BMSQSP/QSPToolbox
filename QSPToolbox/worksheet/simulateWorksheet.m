function myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions)
% This function "crosses" the VPs in a worksheet plus alterations from 
% applied axes with the interventions
% defined int the worksheet and the runs each simulation
% Simulations are run in parallel, if PCT is available.
%
% ARGUMENTS
% myWorksheet:       A worksheet to simulate
% mySimulateOptions: an instance of a simulateOptions object.  If not
%                    provided, simulateWorksheet will use default values.
% 
% RETURNS
% myWorksheet:       the original worksheet plus
%                    the simulation results
%

flagContinue = true;

% First check input arguments
if nargin > 2
    warning(['Too many input arguments to',mfilename,'. Require: myWorksheet; Optional: simulateOptions.'])
    flagContinue = false;
elseif nargin > 1
    flagContinue = true;     
elseif nargin > 0
    mySimulateOptions = simulateOptions;   
    flagContinue = true; 
else
    warning(['Insufficient input arguments to ',mfilename,'. Requires at least: myWorksheet; Optional: simulateOptions.'])
    flagContinue = false;    
end

if flagContinue
    passCheck = mySimulateOptions.verify(myWorksheet);
    if ~passCheck
        warning(['Specified simulation options for ',mfilename,' not valid.'])
        flagContinue = false;        
    end
end

if flagContinue
    % We will assume compiled models are still valid.
    if ~(isequal(class(myWorksheet.compiled.model),'SimBiology.export.Model'))
        warning(['No exported model associated with myWorksheet prior to ',mfilename,', exporting and accelerating.'])
        myWorksheet = compileModel(myWorksheet,true);
        if ~(isequal(class(myWorksheet.compiled.model),'SimBiology.export.Model'))
             warning(['Unable to compile (export & accelerate) model associated with myWorksheet in ',mfilename,'.'])
             flagContinue = false;
        end
    elseif ~myWorksheet.compiled.model.isAccelerated
        warning(['No accelerated model associated with myWorksheet prior to ',mfilename,', accelerating.'])        
        myWorksheet = compileModel(myWorksheet,true);
        if ~(isequal(class(myWorksheet.compiled.model),'SimBiology.export.Model'))
             warning(['Unable to compile (export & accelerate) model associated with myWorksheet in ',mfilename,'.'])
             flagContinue = false;
        end        
    end
end


% Now verify input worksheet
if flagContinue
    vpIDs = getVPIDs(myWorksheet);
    interventionIDs = getInterventionIDs(myWorksheet);
    nVPs = length(vpIDs);
    nInterventions = length(interventionIDs);
    if ((nVPs < 1) || (nInterventions < 1))
        warning(['Insufficient VPs and interventions for ',mfilename,'.'])
        flagContinue = false;
    end
    % At least we can check the specified variables are
    % recongized as elements.  We won't make 'time' requires as specified 
    % but we will add this at the end in regardless.
    if(sum(ismember(myWorksheet.simProps.saveElementResultIDs,myWorksheet.compiled.elements(:,1))) < length(myWorksheet.simProps.saveElementResultIDs))
        warning(['Unable to identify all indicated saveElementResultIDs as elements in myWorksheet in call to ',mfilename,'.'])
        flagContinue = false;
    end    
    
end
    
    
if flagContinue
    
    % We reset the random number generator if specified
    % this will be relevant if we are generating initial guesses
    % for optimization
    if mySimulateOptions.intSeed > -1
        rng(mySimulateOptions.intSeed, 'twister');
    end            
    
    filterFailedRunVPs = mySimulateOptions.filterFailedRunVPs;
    rerunExisting = mySimulateOptions.rerunExisting; 
    optimizeAxisIDs = mySimulateOptions.optimizeAxisIDs;
    optimizeType = mySimulateOptions.optimizeType;
    optimizeMaxIter = mySimulateOptions.optimizeMaxIter;
    
    % We could force re-run for optimization.
    %if sum(ismember({'fmincon'},optimizeType)) > 0
    %    rerunExisting = true;
    %else
    %    rerunExisting = mySimulateOptions.rerunExisting;
    %end
    responseTypeID = mySimulateOptions.responseTypeID;
    % Can only write slices with the parfor loop counter,
    % so simResults needs to be 1-D rather than nInterventions X nVPs
    simResults = cell(1, nVPs * nInterventions);
    % We also need a vector to ultimately track which simulations we will
    % run. We include an option to avoid rerunning all, which is useful
    % for some calling functions.
    nSimulations = nInterventions * nVPs;
    flagRunSim = ones(1, nSimulations);
    flagRunVP = ones(1, nVPs);
    [nInterventionResults, nVPResults] = size(myWorksheet.results);
    if ((nInterventionResults == nInterventions) && (nVPResults == nVPs) && (sum(ismember({'gacohort','psocohort'},optimizeType))<1))
        if ~(rerunExisting)
            myResultClasses = cellfun(@class,myWorksheet.results, 'UniformOutput', false);
            % Results should be stored in a structure, we assume 
            % if a structure is provided then it is a valid result
        	flagRunVP = sum(strcmp(myResultClasses,'struct'),1);
            flagRunVP = ~(flagRunVP == nInterventions);
            flagRunSim = zeros(1, nSimulations);
            for vpCounter = 1 : nVPs
                % increment intervention first, then VP
                if flagRunVP(vpCounter)
                    simulationStartIndex = 1 + (vpCounter - 1) * nInterventions;
                    simulationEndIndex = nInterventions + (vpCounter - 1) * nInterventions;                    
                    flagRunSim(simulationStartIndex:simulationEndIndex) = ones(1,nInterventions);
                %else
                %    flagRunSim = cat(2, flagRunSim, zeros(1, nInterventions));
                end
            end
        end
        % This may consume some nontrivial memory with larger
        % runs, so clear it out.
        clear myResultClasses
    else
        % Otherwise, we don't know which result is paired with which
        % VP & intervention.  We should clear everything out and force
        % a re-run
        myWorksheet.results = cell(nInterventions, nVPs);
    end

    % We need to get the parameters ready before 
    % calling parfor
    [updateValues, updateDoses] = getSimulateValuesDoses(myWorksheet, flagRunVP, flagRunSim);
 
    % In case we will call an optimization, we will need to also
    % get the variable/element indices and bounds along the axes
    if sum(ismember({'fmincon','ga','pso','gacohort','psocohort'},optimizeType)) > 0 
        allAxisIDs = getAxisDefIDs(myWorksheet);
        nOptimizeAxes = length(optimizeAxisIDs);
        % indicesForVaried will be N axes X M interventions
        % as interventions may
        % may overwrite some of the parameters included in the axes, which
        % we will need to account for in the optimization process        
        indicesForVaried = cell(nOptimizeAxes, nInterventions);
        % We make an additional cell for bounds 
        boundsForVaried = cell(nOptimizeAxes, nInterventions);
        axisScale = cell(nOptimizeAxes, 1);
        for axisCounter = 1 : nOptimizeAxes
            myAxisDefID = optimizeAxisIDs{axisCounter};
            curAxisDef = getAxisDef(myWorksheet, myAxisDefID);
            axisIndices = nan(1,0);
            axisBounds = cell(1,0);
            axisScale{axisCounter} = curAxisDef.scale;
            for eCounter = 1 : length(curAxisDef.elementNames)
                curElementName = curAxisDef.elementNames{eCounter};
                curElementType = curAxisDef.elementTypes{eCounter};
                axisBounds = [axisBounds,curAxisDef.bounds{eCounter}];
                theIndices = find(ismember(myWorksheet.compiled.elements(:,1), curElementName) & ismember(myWorksheet.compiled.elements(:,2), curElementType));
                if length(theIndices) > 0
                    axisIndices = cat(2,axisIndices,theIndices);
                    % It should not be possible not to find the axis index
                end
            end    
            for interventionCounter = 1 : nInterventions
                curIntervention = myWorksheet.interventions{interventionCounter};
                [nrows, ncols] = size(curIntervention);
                interventionVariants = extractInterventionTypeElements(curIntervention, 'VARIANT');
                interventionElementValues = flattenVariantstoElements(myWorksheet, interventionVariants);
                checkedAxisIndices = nan(1,0);
                checkedAxisBounds = cell(1,0);
                for eCounter = 1 : length(axisIndices)
                    indexToCheck = axisIndices(eCounter);
                    nameToCheck = myWorksheet.compiled.elements{indexToCheck, 1};
                    typeToCheck = myWorksheet.compiled.elements{indexToCheck, 2};
                    nOverwriteIndices = sum(ismember(interventionElementValues(:,1),nameToCheck) & ismember(interventionElementValues(:,2),typeToCheck));
                    if nOverwriteIndices < 1
                        checkedAxisIndices = cat(2, checkedAxisIndices, indexToCheck);
                        checkedAxisBounds = [checkedAxisBounds, axisBounds{eCounter}];
                    % else
                        % If we don't pass the check, the intervention will
                        % overwrite the axis element, so we don't add this
                        % element as a variable
                    %     1;
                    end
                end
                indicesForVaried{axisCounter,interventionCounter} = checkedAxisIndices;
                boundsForVaried{axisCounter,interventionCounter} = checkedAxisBounds;
            end
        end
        optAxesCoefs = cell(1,nVPs);
    end
    
    % Specify simulation length and sample times
    exportedModel = myWorksheet.compiled.model;
    exportedModel.SimulationOptions.OutputTimes = myWorksheet.simProps.sampleTimes;
    exportedModel.SimulationOptions.StopTime = max(myWorksheet.simProps.sampleTimes);
    exportedModel.SimulationOptions.AbsoluteTolerance = myWorksheet.simProps.absoluteTolerance;
    exportedModel.SimulationOptions.RelativeTolerance = myWorksheet.simProps.relativeTolerance;    
    exportedModel.SimulationOptions.AbsoluteToleranceScaling = myWorksheet.simProps.absoluteToleranceScaling;    
    exportedModel.SimulationOptions.AbsoluteToleranceStepSize = myWorksheet.simProps.absoluteToleranceStepSize;    
    exportedModel.SimulationOptions.MaximumWallClock = myWorksheet.simProps.maximumWallClock;
    exportedModel.SimulationOptions.SolverType = myWorksheet.simProps.solverType;
    % It doesn't appear the exportedModel has a property similar to states 
    % to log in the configset object, so this will be done manually.
    mySaveElementResultIDs = myWorksheet.simProps.saveElementResultIDs;
    
    % As a precaution, restart any existing parallel 
    % pools 
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end
    parpool;
	% in 2017a this message can cause issues, especially when
	% trying some of the global optimization methods
	% You may want this on for some applications,
	% but it is disabled here.
	parfevalOnAll(gcp(), @warning, 0, 'off', 'SimBiology:Simulation:NonFiniteRHS');

    simResults = cell(1,nSimulations);
    if sum(ismember({'none'},optimizeType)) > 0
        parfor simulationCounter = 1 : nSimulations
        % for simulationCounter = 1 : nSimulations    
            if flagRunSim(simulationCounter) > 0
                % increment intervention first, then VP
                % so we can complete simulations for a given VP
                vpCounter = floor((simulationCounter-1)/(nInterventions)) + 1;
                interventionCounter = simulationCounter - (vpCounter - 1) * nInterventions;
                % Retrieve the pre-specified parameter alterations & doses
                % for the current run
                currentModelValues = updateValues(interventionCounter, vpCounter,:);
                currentModelValues = reshape(currentModelValues,[],1);
                curDoses = updateDoses{interventionCounter};
                % Simulate with specified sampling interval
                try
                    simData = simulate(exportedModel, currentModelValues, curDoses);
                    simData = convertSimData(simData);
                    % Reformat results into data table format
                    % TODO: we will want to add some more checks to the saveElementResultIDs
                    if length(mySaveElementResultIDs) > 0
                        keepIndices = find(ismember(simData.Names,['time',mySaveElementResultIDs]));
                        simData.Names = simData.Names(keepIndices);
                        simData.Data = simData.Data(:,keepIndices);
                    end
                    flagSimData = verifySimData(simData, exportedModel.SimulationOptions.StopTime);
                    % If the simulation doesn't complete the entire
                    % specified length of simulated time, then we need to
                    % flag this.  For now, we only keep results
                    if ~flagSimData
                        simResults{simulationCounter} = [];
                    else
                        simResults{simulationCounter} = simData;
                    end
                catch
                    % If the simulation fails, we need to move on
                    simResults{simulationCounter} = [];
                end
            else
                simResults{simulationCounter} = [];
            end
        end
        simResults = reshape(simResults, [nInterventions, nVPs]); 
        if sum(flagRunVP) < nVPs
            myWorksheet.results(:,find(flagRunVP)) = simResults(:,find(flagRunVP));
        else
            myWorksheet.results = simResults;
        end 
    elseif sum(ismember({'fmincon','ga','pso'},optimizeType)) > 0
    % elseif sum(ismember({'fmincon','ga'},optimizeType)) > 0
    	mySeedCohort = getVPCoeffs(myWorksheet);
        for vpCounter = 1 : nVPs
            coeffsToVary = nan(nOptimizeAxes,1);
            optimizeAxisIndices = find(ismember(allAxisIDs,optimizeAxisIDs));
            if mySimulateOptions.optimizeSeedExisting
                coeffsToVary = transpose(mySeedCohort(optimizeAxisIndices,vpCounter));
                initialPopulation = [coeffsToVary;lhsdesign(mySimulateOptions.optimizePopSize-1,nOptimizeAxes)];
            else
                coeffsToVary = [];
                initialPopulation = lhsdesign(mySimulateOptions.optimizePopSize,nOptimizeAxes);
            end
            if sum(ismember({'ga'},optimizeType)) > 0
                optimOptions = gaoptimset;
                optimOptions.Display = 'diagnose';
                % Time limit in s.  We will probably want to move this
                % to a user-defined option.
                optimOptions.TimeLimit = mySimulateOptions.optimizeTimeLimit;
                optimOptions.PopulationSize = mySimulateOptions.optimizePopSize;
                % We may want to switch this back to uniform,
                % but I was working with one problem that would benefit
                % from the adaptfeasible.
                optimOptions.MutationFcn = {@mutationuniform, 1.5/nOptimizeAxes};
                %optimOptions.MutationFcn = {@mutationadaptfeasible, 1.5/nOptimizeAxes};
                optimOptions.PopInitRange = cat(1,zeros(1,nOptimizeAxes),ones(1,nOptimizeAxes));
                optimOptions.InitialPopulation = initialPopulation;  
                optimOptions.TolFun = 1E-3;
                optimOptions.UseParallel = true;
                if optimizeMaxIter > -1
                    optimOptions.Generations = optimizeMaxIter;
                end
                % I was having issues getting hybridopts-fmincon to run in
                % parallel consistently, so this is disabled for now.
                % hybridopts = optimoptions('fmincon');
                % hybridopts.Display = 'iter';
                % hybridopts.UseParallel = true;
                % hybridopts.Algorithm = 'active-set';
                % optimOptions.HybridFcn = {@fmincon,hybridopts};                
            elseif sum(ismember({'pso'},optimizeType)) > 0
                optimOptions = optimoptions('particleswarm');
                optimOptions.Display = 'iter';
                % Time limit in s.  We will probably want to move this
                % to a user-defined option.
                optimOptions.MaxTime = mySimulateOptions.optimizeTimeLimit;
                optimOptions.SwarmSize = mySimulateOptions.optimizePopSize;
                optimOptions.TolFun = 1E-3;
                optimOptions.UseParallel = true;                
                optimOptions.InitialSwarm = initialPopulation;
                % Since the bounds are finite, the InitialSwarmSpan option 
                % is ignored by the particle swarm algorithm so we don't
                % need to worry about it.
                if optimizeMaxIter > -1
                    optimOptions.MaxIterations = optimizeMaxIter;
                end                    
            else 
                optimOptions = optimoptions('fmincon');
                optimOptions.Display = 'iter';
                optimOptions.UseParallel = true;
                % Not an option for fmincon
                % optimOptions.TimeLimit = mySimulateOptions.optimizeTimeLimit;
                % Algorithm options that do not require supplying a gradient function:
                % 'interior-point' (default, large scale)
                % 'sqp'
                % 'active-set'
                % optimOptions.Algorithm = 'sqp';
                % optimOptions.Algorithm = 'active-set';
                % Testing on 12/25/15 with the ADC platform suggested they
                % perform comparably, perhaps with activeset behaving the
                % best by a small margin.
                optimOptions.Algorithm = 'active-set';
                if optimizeMaxIter > -1
                    optimOptions.MaxIterations = optimizeMaxIter;
                end                    
            end            
            myVPID = vpIDs{vpCounter};
            reducedWorksheet = copyWorksheet(myWorksheet,{myVPID}, false);
            if flagRunVP(vpCounter) > 0
                vpElementDefaultValues = updateValues(:, vpCounter,:);
                [nElements, dummy] = size(myWorksheet.compiled.elements);
                elementValuesAcrossInterventions = nan(nElements, nInterventions);
                for interventionCounter = 1 : nInterventions
                    interventionElementDefaultValues = vpElementDefaultValues(interventionCounter,1,:);
                    elementValuesAcrossInterventions(:,interventionCounter) = reshape(interventionElementDefaultValues,[],1);
                end
                % I ran into issues with optimization in the parfor loop,
                % SimBiology stopped recognizing these as belonging to the
                % original exported model.  This issue did not exist in a for
                % loop.  So this results in some duplication of coding, and
                % a little extra computation, but I've re-exported the doses
                % here.
                theDoses = updateDoses(:);
                for interventionCounter = 1 : nInterventions
                    curDoses = theDoses(interventionCounter);
                    if ~isempty(curDoses{1})
                        for theDoseCounter = 1 : length(curDoses{1})
                            curDose = curDoses{1}(theDoseCounter);
                            if ~isempty(curDose)
                                curDoses{1}(theDoseCounter) = getdose(exportedModel,curDose.Name);
                            else
                                curDoses{1}(theDoseCounter) = curDose;
                            end
                        end
                    end
                    theDoses(interventionCounter)=curDoses;
                end
                % We may want to add a cehck here in case
                anonymousFunction = @(x)runModelforOptimization(exportedModel, reducedWorksheet, elementValuesAcrossInterventions, theDoses, x, indicesForVaried, boundsForVaried, axisScale, responseTypeID);
                if sum(ismember({'ga'},optimizeType)) > 0
                    [optResult,fVal,exitFlag,output] = ga(anonymousFunction,nOptimizeAxes,[],[],[],[],zeros(nOptimizeAxes,1),ones(nOptimizeAxes,1),[],optimOptions);
                    optResult = transpose(optResult);
                elseif sum(ismember({'pso'},optimizeType)) > 0
                    [optResult,fVal,exitFlag,output] = particleswarm(anonymousFunction,nOptimizeAxes,zeros(nOptimizeAxes,1),ones(nOptimizeAxes,1),optimOptions);
                    optResult = transpose(optResult);                    
                else
                    [optResult,fVal,exitFlag,output] = fmincon(anonymousFunction,coeffsToVary,[],[],[],[],zeros(1,nOptimizeAxes),ones(1,nOptimizeAxes),[],optimOptions);
                end
                % disp(exitFlag)
                % Only exitFlag -1, -2 are really bad according to
                % documentation for fmincon and ga.  1 is best.  We'll be relatively
                % loose in accepting outcomes 0, 1, 2...
                if ((exitFlag ~= -1) && (exitFlag ~= -2))
                    optAxesCoefs{vpCounter} = optResult;
                    
                else
                    optAxesCoefs{vpCounter} = [];
                end
            else
                optAxesCoefs{vpCounter} = [];
            end
        end
        for vpCounter = 1 : nVPs
            if length(optAxesCoefs{vpCounter}) > 0
                if flagRunVP(vpCounter) > 0
                    for axisCounter = 1 : nOptimizeAxes
                        myAxisDefID = optimizeAxisIDs{axisCounter};
                        axisIndex = find(ismember(allAxisIDs, myAxisDefID));
                        myWorksheet.axisProps.axisVP.coefficients(axisIndex,vpCounter) = optAxesCoefs{vpCounter}(axisCounter);
                    end
                end
            else
                % If we couldn't optimize, and a result already exists,
                % we will keep the original axis coefficients for the
                % VP and the original results.
                % We can set flagRunVP(vpCounter) here to
                % zero, with the knowledge that since this was
                % an optimization run the simulation will be re-run
                % if results are missing anyway.
                if ~(rerunExisting)
                    flagRunVP(vpCounter) = 0;
                else
                    flagRunVP(vpCounter) = 1;
                end
            end
        end
    elseif sum(ismember({'gacohort','psocohort'},optimizeType)) > 0
        optimizeAxisIndices = find(ismember(allAxisIDs,optimizeAxisIDs));
        if mySimulateOptions.optimizeSeedExisting
            mySeedCohort = getVPCoeffs(myWorksheet);
            mySeedCohort = transpose(mySeedCohort(optimizeAxisIndices,:));
            if nVPs < mySimulateOptions.optimizePopSize
                mySeedCohort = [mySeedCohort;lhsdesign(mySimulateOptions.optimizePopSize-nVPs,nOptimizeAxes)];
            elseif nVPs > mySimulateOptions.optimizePopSize
                mySeedCohort = datasample(mySeedCohort,mySimulateOptions.optimizePopSize,1);
            end
        else
            mySeedCohort = lhsdesign(mySimulateOptions.optimizePopSize,nOptimizeAxes);
        end
                
        if sum(ismember({'gacohort'},optimizeType)) > 0
            optimOptions = gaoptimset;
            %optimOptions.HybridFcn = {@fmincon};
            % We weill leave this on a verbose setting for now,
            % but we may want to change later to 'off'
            optimOptions.Display = 'diagnose';
            optimOptions.TolFun = 1E-3;
            % Default optimOptions.PopulationSize: '50 when numberOfVariables
            % <= 5, else 200' 
            optimOptions.PopulationSize = mySimulateOptions.optimizePopSize;
            optimOptions.EliteCount = round(0.025*optimOptions.PopulationSize);
            % We may want to reduce the CrossoverFraction to introduce more
            % mutations.  Defauly is 0.8.
            optimOptions.CrossoverFraction = 0.6;
            % Default: @fitscalingrank
            % This did not work as well optimOptions.FitnessScalingFcn = @fitscalingprop;
            optimOptions.FitnessScalingFcn = @fitscalingrank;
            % Default optimOptions.Generations: '100*numberOfVariables'
            % For mutations, MATLAB will default to mutationadaptfeasible
            % for consintrained minimization.  However, we'd like to make
            % sure that we are doing a good job exploring the allowed search
            % ranges. So we try to implement mutationunform, which uses  
            % optimOptions.PopInitRange to set the boundaries, the default
            % rate in MATLAB is 0.01.
            optimOptions.MutationFcn = {@mutationuniform, 1.5/nOptimizeAxes};
            optimOptions.PopInitRange = cat(1,zeros(1,nOptimizeAxes),ones(1,nOptimizeAxes));
            %popSizeToReturn = 200;
            optimOptions.InitialPopulation = mySeedCohort;
            optimOptions.UseParallel = true;
            optimOptions.TimeLimit = mySimulateOptions.optimizeTimeLimit;
            if optimizeMaxIter > -1
                optimOptions.Generations = optimizeMaxIter;
            end
        else
            optimOptions = optimoptions('particleswarm');
            optimOptions.Display = 'iter';
            % Time limit in s.  We will probably want to move this
            % to a user-defined option.
            optimOptions.MaxTime = mySimulateOptions.optimizeTimeLimit;
            optimOptions.SwarmSize = mySimulateOptions.optimizePopSize;
            optimOptions.TolFun = 1E-3;
            optimOptions.UseParallel = true;
            optimOptions.InitialSwarm =  mySeedCohort;
            if optimizeMaxIter > -1
                optimOptions.MaxIterations = optimizeMaxIter;
            end            
        end
        % We will specify the seed population, so we can
        % can just use the first worksheet VP as a template for cohort
        % optimization runs
        myVPID = vpIDs{1};
        reducedWorksheet = copyWorksheet(myWorksheet,{myVPID}, false);
        vpElementDefaultValues = updateValues(:, 1,:);
        [nElements, ~] = size(myWorksheet.compiled.elements);
        elementValuesAcrossInterventions = nan(nElements, nInterventions);
        for interventionCounter = 1 : nInterventions
            interventionElementDefaultValues = vpElementDefaultValues(interventionCounter,1,:);
            elementValuesAcrossInterventions(:,interventionCounter) = reshape(interventionElementDefaultValues,[],1);
        end
        theDoses = updateDoses(:);
        for interventionCounter = 1 : nInterventions
            curDoses = theDoses(interventionCounter);
            if ~isempty(curDoses{1})
                for theDoseCounter = 1 : length(curDoses{1})
                    curDose = curDoses{1}(theDoseCounter);
                    if ~isempty(curDose)
                        curDoses{1}(theDoseCounter) = getdose(exportedModel,curDose.Name);
                    else
                        curDoses{1}(theDoseCounter) = curDose;
                    end
                end
            end
            theDoses(interventionCounter)=curDoses;
        end
        anonymousFunction = @(x)runModelforOptimization(exportedModel, reducedWorksheet, elementValuesAcrossInterventions, theDoses, x, indicesForVaried, boundsForVaried, axisScale, responseTypeID);
        if sum(ismember({'gacohort'},optimizeType)) > 0
            % Note that for 'gacohort', MATLAB returns:
            % "Warning: You are using 'mutationuniform' mutation function for
            % constrained minimization.  Solution may be infeasible; use
            % '@mutationadaptfeasible' function for constrained minimization." 
            % However, we have strictly constant bounds around the edges
            % so the mutationuniform function should note generate infeasible
            % points.  So we disable this warning.   
            warning('off','globaloptim:constrvalidate:unconstrainedMutationFcn');
            [optResult,fVal,exitFlag,output,finalPopulation, scores] = ga(anonymousFunction,nOptimizeAxes,[],[],[],[],zeros(nOptimizeAxes,1),ones(nOptimizeAxes,1),[],optimOptions);
            warning('on','globaloptim:constrvalidate:unconstrainedMutationFcn');
        % 'psocohort' is the other case.
        else
            % The custom tweaked version is needed because MathWork's
            % version doesn't provide the swarm.
            [optResult,fVal,exitFlag,finalPopulation] = particleSwarmCohortWrapper(anonymousFunction,nOptimizeAxes,zeros(nOptimizeAxes,1),ones(nOptimizeAxes,1),optimOptions);
        end
        finalPopulation = transpose(finalPopulation);
        [~,finalPopSize] = size(finalPopulation);
        allAxisDefIDs = getAxisDefIDs(myWorksheet);
        firstVP = getVP(reducedWorksheet, myVPID);
        namedCounter = 0;
        newVPIDs = cell(1,finalPopSize);
        newVPIDs(1,:)={''};
        newVPVariants = cell(1,finalPopSize);
        while namedCounter < finalPopSize
            testNumber = 1;
            namedCounter = namedCounter + 1;
            nameTry = [myVPID,'_',num2str(testNumber)];
            while ((sum(ismember(vpIDs,nameTry)) + sum(ismember(newVPIDs,nameTry))) > 0)
                testNumber = testNumber + 1;
                nameTry = [myVPID,'_',num2str(testNumber)];
            end
            newVPIDs{1,namedCounter} = nameTry;
            newVPVariants{1,namedCounter} = firstVP.variants;
            % newVPAxes{1,namedCounter} = reducedWorksheet.axisProps.axisVP(:,1);
        end
        % We also need to organize the VP coefficient matrix
        newVPCoeffs = getVPCoeffs(reducedWorksheet);
        newVPCoeffs = newVPCoeffs*ones(1,finalPopSize); 
        for axisCounter = 1 : length(allAxisDefIDs)            
            axisID = allAxisDefIDs(axisCounter);
            if sum(ismember(optimizeAxisIDs,axisID)) > 0
                optimizeAxisIndex = find(ismember(optimizeAxisIDs, axisID));
                newVPCoeffs(axisCounter,:) = finalPopulation(optimizeAxisIndex,:);
            end
        end    
        reducedWorksheet = createVPs(reducedWorksheet, newVPIDs, newVPVariants, newVPCoeffs);
        reducedWorksheet = removeVPs(reducedWorksheet,{myVPID});
        myWorksheet = reducedWorksheet;
        nVPs = finalPopSize;
        nSimulations = finalPopSize * nInterventions;
        simResults = cell(1, nSimulations);
        flagRunSim = ones(1, nSimulations);
        flagRunVP = ones(1, finalPopSize);        
    end
    
    simResults = reshape(simResults, [nInterventions, nVPs]);
    % Clean up the worker pool
    delete(gcp)
    if sum(flagRunVP) < nVPs
        myWorksheet.results(:,find(flagRunVP)) = simResults(:,find(flagRunVP));
    else
        myWorksheet.results = simResults;
    end
    
    
    if sum(ismember({'fmincon','ga','pso','gacohort','psocohort'},optimizeType)) > 0
        % If we just ran an optimization, update the results stored in the
        % worksheet with the new VPs
        newSimulateOptions = mySimulateOptions;
        % We will filter failed VPs in this current run, not the additional
        % call.
        newSimulateOptions.optimizeType = 'none';
        newSimulateOptions.responseTypeID = ''; 
        newSimulateOptions.optimizeAxisIDs = {};
        newSimulateOptions.filterFailedRunVPs = false;
        % saveElementResultIDs will be passed along with the worksheet
        myWorksheet = simulateWorksheet(myWorksheet, newSimulateOptions);
    end
 
    % Filter VPs with failed results?  These are potentially useful as a
    % diagnostic
    if filterFailedRunVPs
        filterVPids = cell(1,0);
        allVPIDs = getVPIDs(myWorksheet);
        for vpCounter = 1 : nVPs
            if sum(cellfun('isempty',myWorksheet.results(:,vpCounter))) > 0
                filterVPids{length(filterVPids)+1} = allVPIDs{vpCounter};
            end
        end
        myWorksheet = removeVPs(myWorksheet, filterVPids);
    end

else
    warning(['Could not complete ',mfilename,'. Returning original worksheet.'])
end
end