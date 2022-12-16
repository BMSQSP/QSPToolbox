function [myWorksheet, mySimulateTolerances, simStats] = simulateWorksheetIterateTolerance(myWorksheet, mySimulateOptions, opts)
    % This function will iterate with decreased tolerance and try
    % to force all of the VPs in a worksheet to complete.    
    %
    % USAGE
    %   [myWorksheet, mySimulateTolerances, simStats] = simulateWorksheetIterateTolerance(myWorksheet, mySimulateOptions, <optional keyword arguments...>)   
    %
    % ARGUMENTS
    %
    %   myWorksheet:
    %       A worksheet to simulate.
    %
    %   mySimulateOptions: optional, simulateOptions object
    %       If not provided, this function will use default values.
    %
    %   passThreshold: keyword, optional
    %       Percentage of successful simulations to stop the iterations.
    %       Default: 1.00.
    %
    %   maxAbsTol: keyword, optional
    %       Maximum absolute tolerance to try. Default: 1e-4.
    %
    %   minAbsTol: keyword, optional
    %       Minimum absolute tolerance to try. Default: 1e-19.
    %
    %   maxRelTol: keyword, optional
    %       Maximum relative tolerance to try. Default: 1e-4.
    %
    %   minRelTol: keyword, optional
    %       Minimum relative tolerance to try. Default: 1e-15.
    %
    %   absTolStepSize: keyword, optional
    %       Absolute Tolerance Step Size for use with Absolute Tolerance
    %       Scaling. Default: 1e-10.
    %
    %   numAbsTolLevels: keyword, optional
    %       Number of grid points per absolute or relative tolerance axis. The
    %       total number of (rel, abs) combinations will be
    %       numAbsTolLevels^2. Default: 5.
    %
    %   tryHigherAbsTols: keyword, optional
    %       Whether to attempt absolute tolerance levels that are higher
    %       than what is in myWorksheet.simProps. Default: false.
    %
    %   printStatus: keyword, optional
    %       Logical value for whether to display intermediate status. Default:
    %       false.
    %
    % RETURNS
    %
    %  myWorksheet:           the original worksheet updated with the
    %                           the simulation results
    %
    %  mySimulateTolerances:  the simProps settings from the worksheet
    %                           (optional, struct)
    %
    %  simStats:              table of simulation success rate per
    %                         combination of tolerances attempted.
    %
    %  
    %
    % SEE ALSO
    %   simulateWorksheet, simulateOptions
    arguments
        myWorksheet
        mySimulateOptions (1,1) simulateOptions = simulateOptions()
        opts.passThreshold (1,1) double = 1.00;
        opts.maxAbsTol (1,1) double = 1.0E-6; % 1e-6 is the default value for AbsoluteTolerance in SimBiology
        opts.minAbsTol (1,1) double = 1.0E-19;
        opts.maxRelTol (1,1) double = 1.0E-4;
        opts.minRelTol (1,1) double = 1.0E-15;
        opts.absTolStepSize (1,1) double = 1.0E-10;
        opts.numRelTolLevels (1,1) {mustBeInteger, mustBePositive} = 1;
        opts.numAbsTolLevels (1,1) {mustBeInteger, mustBePositive} = 5;
        opts.printStatus (1,1) logical = false;
        opts.tryHigherAbsTols (1,1) logical = false;
    end

    resimulateInitial = mySimulateOptions.rerunExisting;

    flagContinue = true;
    if flagContinue
        passCheck = mySimulateOptions.verify(myWorksheet);
        if ~passCheck
            warning(['Specified simulation options for ',mfilename,' not valid.'])
            flagContinue = false;
        end
        if mySimulateOptions.filterFailedRunVPs
            mySimulateOptions.filterFailedRunVPs = false;
            flagFilterAtEnd = true;
        else
            flagFilterAtEnd = false;
        end
    end

    % Just initialize a cell array to return if requested, will
    % update later
    mySimulateTolerances = cell(1, 1);

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
        mySimulateTolerances = cell(nInterventions, nVPs);
    end

    if ~flagContinue
        error(['Cannot finish ', mfilename()]);
    end

    if resimulateInitial
        myWorksheet.results = {};
    end

    % As a precaution, restart any existing parallel
    % pools
    if mySimulateOptions.poolRestart
        if ~isempty(gcp('nocreate'))
            delete(gcp);
        end
    end
    if isempty(gcp('nocreate'))
        % First check the default number of workers, if needed
        mySimulateOptions = checkNWorkers(mySimulateOptions);
        parpool(mySimulateOptions.clusterID,mySimulateOptions.nWorkers,'SpmdEnabled',false);
    end

    % Stash the worksheet's original tolerance
    originalSimProps = myWorksheet.simProps;

    % Set up options for simulateWorksheet
    mySimulateOptions.poolRestart = false;
    mySimulateOptions.poolClose = false;
    mySimulateOptions.rerunExisting = false;
    mySimulateOptions.filterFailedRunVPs = false;

    % Generate the tolerance levels to try, these
    % tolerances are uniform in logspace
    relTolLevels = logspace( ...
        log10(opts.minRelTol), ...
        log10(opts.maxRelTol), ...
        opts.numRelTolLevels ...
        );
    % We also add the tolerances specified in the input worksheet
    relTolLevels = unique([relTolLevels, myWorksheet.simProps.relativeTolerance]);
    relTolLevels = sort(relTolLevels, 'descend');

    % For absolute tolerance, we first attempt simulation with reduced levels,
    % then with increased levels for the remaining unsuccessful VPs
    lowerAbsTolLevels = logspace( ...
            log10(opts.minAbsTol), ...
            log10( ...
                    max( ...
                    opts.minAbsTol, ...
                    min(myWorksheet.simProps.absoluteTolerance, opts.maxAbsTol) ...
                    ) ...
            ),...
            opts.numAbsTolLevels...
        );
    lowerAbsTolLevels = sort(lowerAbsTolLevels, 'descend');
    higherAbsTolLevels = logspace( ...
            log10( ...
                    min( ...
                    opts.maxAbsTol, ...
                    max(myWorksheet.simProps.absoluteTolerance, opts.minAbsTol)...
                    )...
                ), ...
            log10(opts.maxAbsTol),...
            opts.numAbsTolLevels...
        );
    higherAbsTolLevels = sort(higherAbsTolLevels, 'ascend');    

    % Generate all combinations of relative tolerance and absolute
    % tolerance levels, ordered by quadrants

    % First quadrant: reducing both relative Tol and abs Tol
    [absTolLevels, relTolLevels] = ndgrid(lowerAbsTolLevels, relTolLevels);
    firstQuadrantTolOptions = [reshape(relTolLevels, [], 1) reshape(absTolLevels, [], 1)];

    % Second quadrant: reducing relative Tol but try increasing abs Tol
    % This quadrant is only explored if user set keyword argument
    % tryHigherAbsTols to true
    if opts.tryHigherAbsTols
        [absTolLevels, relTolLevels] = ndgrid(higherAbsTolLevels, relTolLevels);
        secondQuadrantTolOptions = [reshape(relTolLevels, [], 1) reshape(absTolLevels, [], 1)];
    else
        secondQuadrantTolOptions = zeros(0, 2);
    end    

    tolOptions = [
        firstQuadrantTolOptions;
        secondQuadrantTolOptions];
    tolOptions = unique(tolOptions, 'stable');

    % Set up solver options
    if strcmp(myWorksheet.simProps.solverType, 'sundials')
        solverOptions = ["sundials", "ode15s"];
    else
        solverOptions = ["ode15s", "sundials"];
    end

    % Set up table to store performance metrics
    simStats = table( ...
        'Size',...
        [0, 8],...
        'VariableNames', ...
        {'absTol', 'relTol', 'absScaling', 'absTolStepSize', 'solver', 'numSimulations', 'numFailures', 'wallclockTime'}, ...
        'VariableTypes', ...
        {'double', 'double', 'logical', 'double', 'string', 'int32', 'int32', 'double'}...
        );

    totalVPCount = length(myWorksheet.vpDef);    
    for solver = solverOptions
        for tolOptionIter = 1:length(tolOptions)
            relTol = tolOptions(tolOptionIter, 1);
            absTol = tolOptions(tolOptionIter, 2);

            tic
            [status, nSimFail] = simulateWithNewTolAndSolver(relTol, absTol, solver);
            wcTime = toc;

            if size(simStats, 1) > 0
                simStats{end, 'numFailures'} = nSimFail;
            end

            % Add new row to simulation performance table
            newRow = struct('absTol', absTol,...
                'relTol', relTol,...
                'absScaling', true, ...
                'absTolStepSize', opts.absTolStepSize,...
                'solver', solver,...
                'numSimulations', nSimFail,...
                'wallclockTime', wcTime);
            simStats = [simStats; struct2table(newRow)];            

            if status == 1
                break;
            end
        end
    end

    % Restore the original simulation tolerances
    myWorksheet.simProps = originalSimProps;

    % Clean up the pool, if needed
    if mySimulateOptions.poolClose
        if ~isempty(gcp('nocreate'))
            delete(gcp);
        end
    end

    % Follow desired behavior for failed simulations.
    if nSimFail > 0
        if ~flagFilterAtEnd
            disp(['Not all simulations completed in ',mfilename,'.  ' ...
                'There are still ',num2str(nSimFail),' unsuccessful simulations.  ' ...
                'Returning worksheet with available results and all VPs.'])
        else
            disp(['Not all simulations completed in ',mfilename,'.  ' ...
                'There are still ',num2str(nSimFail),' unsuccessful simulations.  ' ...
                'Returning worksheet with available results and successful VPs.'])
            allVPIDs = getVPIDs(myWorksheet);
            myResultClasses = cellfun(@class,myWorksheet.results, 'UniformOutput', false);
            failFlags = ~(strcmp(myResultClasses,'struct'));
            failVPLogical = sum(failFlags,1) > 0;
            filterVPids = allVPIDs(failVPLogical);
            myWorksheet = removeVPs(myWorksheet, filterVPids);
            mySimulateTolerances = mySimulateTolerances(:,~failVPLogical);
        end
    end

    function [status, nSimFail] = simulateWithNewTolAndSolver(relTol, absTol, solver)
        status = 0;
        solverType = char(solver);
        % Get status of existing prior results and initial tolerances
        failureFlags = findFailedSimulations(myWorksheet);
        nSimFail = sum(sum(failureFlags));

        if opts.printStatus
            disp(['Failure in ', num2str(nSimFail), ...
                ' VPxIntervention combinations.']);
        end

        % We will not attempt to solve again if the fraction of
        % failures is acceptable
        if nSimFail <= (1 - opts.passThreshold)*totalVPCount
            status = 1;
            return;
        end

        % Set tolerance levels
        myWorksheet.simProps.absoluteTolerance = absTol;
        myWorksheet.simProps.relativeTolerance = relTol;
        myWorksheet.simProps.absoluteToleranceScaling = true;
        myWorksheet.simProps.absoluteToleranceStepSize = opts.absTolStepSize;
        myWorksheet.simProps.solverType = solverType;

        % Re-simulate the worksheet
        myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions);

        % Now assign the obtained results to the original
        % worksheet. We only need to re-assign the new simulations,
        % not the already successful ones
        runSimFlags = reshape(failureFlags, 1, []);
        mySimulateTolerances(runSimFlags) = {
            struct( ...
            'sampleTimes', [myWorksheet.simProps.sampleTimes], ...
            'absoluteTolerance', absTol, ...
            'absoluteToleranceStepSize', opts.absTolStepSize, ...
            'relativeTolerance', relTol, ...
            'absoluteToleranceScaling', true, ...
            'maximumWallClock', myWorksheet.simProps.maximumWallClock, ...
            'maxStep', myWorksheet.simProps.maxStep, ...
            'solverType', myWorksheet.simProps.solverType, ...
            'saveElementResultIds', {myWorksheet.simProps.saveElementResultIDs}...
            )
            };
    end
end

function failureFlags = findFailedSimulations(worksheet)
    % Generate a nInterventions x nVPs true/false matrix to flag the
    % simulations that failed
    results = worksheet.results;
    interventionCount = length(worksheet.interventions);
    vpCount = length(worksheet.vpDef);

    failureFlags = true(interventionCount, vpCount);

    if isempty(results)
        return;
    end

    for vpIter = 1:vpCount
        for interventionIter = 1:interventionCount
            if ~isempty(results{interventionIter, vpIter})
                failureFlags(interventionIter, vpIter) = false;
            end
        end
    end
end