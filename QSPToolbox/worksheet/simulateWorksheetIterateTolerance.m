function myWorksheet = simulateWorksheetIterateTolerance(myWorksheet, mySimulateOptions)
% This function will iterate with decreased tolerance or max step and try to
% force all of the VPs in a worksheet to complete.  This function
% will iteratively call simulate worksheet.  
% NOTE 1: This function is provided for convenience to 
%         try to force simulation of all of
%         the provided VPs/parameter vectors to complete, but they 
%         will not necessarily all be completed with the
%         same settings.
% NOTE 2: By default we will
% 		  resimulate the initial VPs unless a mySimulateOptions is 
% 		  provided specifying otherwise.
%
% ARGUMENTS
%  myWorksheet:       A worksheet to simulate.
%  mySimulateOptions: An instance of a simulateOptions object.  If not
%                     provided, this function will use default values.
% 
% RETURNS
% myWorksheet:       the original worksheet plus
%                    the simulation results
%

flagContinue = true;

% First check input arguments
resimulateInitial = true;
if nargin > 2
    warning(['Too many input arguments to',mfilename,'. Require: myWorksheet; Optional: simulateOptions.'])
    flagContinue = false;
elseif nargin > 1
    flagContinue = true;  
	resimulateInitial = mySimulateOptions.filterFailedRunVPs;
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
	if mySimulateOptions.filterFailedRunVPs
        warning(['The filterFailedRunVPs flag for mySimulateOptions in ',mfilename,' should be set to false.  Resetting.'])
        mySimulateOptions.filterFailedRunVPs = false;        
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
	if resimulateInitial
		myWorksheet.results = {};
	end
	mySimulateOptions.rerunExisting = false;
	myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions);
	myResultClasses = cellfun(@class,myWorksheet.results, 'UniformOutput', false);
	nSimFail = sum(sum(~(strcmp(myResultClasses,'struct'))));
	originalAbsTol = myWorksheet.simProps.absoluteTolerance;
	originalRelTol = myWorksheet.simProps.relativeTolerance;
	nRetries = 1;
	resimulateFlag = true;
	% This is clearly empirical, but try 3 resubmissions
	% with decreasing absolute tolerance
	while ((nSimFail > 0) && (resimulateFlag) && (nRetries<=3))
		myWorksheet.simProps.absoluteTolerance = originalAbsTol/(10*nRetries);
		myWorksheet.simProps.relativeTolerance = originalRelTol;
		myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions);
		myResultClasses = cellfun(@class,myWorksheet.results, 'UniformOutput', false);
		% Results should be stored in a structure, we assume 
		% if a structure is provided then it is a valid result
		nSimFail = sum(sum(~(strcmp(myResultClasses,'struct'))));	
		nRetries = nRetries + 1;
		if originalAbsTol/(10*nRetries) < 1E-100
			resimulateFlag = false
		end
	end
	myWorksheet.simProps.absoluteTolerance = originalAbsTol;
	myWorksheet.simProps.relativeTolerance = originalRelTol;
	originalMaxStep = myWorksheet.simProps.maxStep;	
	if isnumeric(originalMaxStep)
		testStepBasis = originalMaxStep;
	else
		testStepBasis = min(diff(myWorksheet.simProps.sampleTimes));
	end
	% If that fails, try adjusting the relative step size
	nRetries = 1;
	resimulateFlag = true;
	while ((nSimFail > 0) && (resimulateFlag) && (nRetries<=3))
		myWorksheet.simProps.maxStep = testStepBasis/(10*nRetries);
		myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions);
		myResultClasses = cellfun(@class,myWorksheet.results, 'UniformOutput', false);
		% Results should be stored in a structure, we assume 
		% if a structure is provided then it is a valid result
		nSimFail = sum(sum(~(strcmp(myResultClasses,'struct'))));	
		nRetries = nRetries + 1;
		if testStepBasis/(10*nRetries) < 1E-9
			resimulateFlag = false;
		end
	end	
	myWorksheet.simProps.maxStep = originalMaxStep;
	if nSimFail > 0
		warning(['Not all simulations completed in ',mfilename,'.  There are still ',num2str(nSimFail),' unsuccessful simulations.  Returning worksheet with current results and all VPs.'])
	end
else
    warning(['Could not complete ',mfilename,'. Returning original worksheet.'])
end
end