% Here, we illustrate workflow for creating VPops using MAPEL
%
% Toolbox initialization
% toolboxPath = 'DRIVE:/DIRs/QSPToolbox';
% addpath(toolboxPath);
initQSPToolbox;

% Load in the worksheet and regenerate results.
myWorksheetFileName = 'example_cohort_worksheet1';
myWorksheet = loadWorksheet(myWorksheetFileName);
mySimulateOptions = simulateOptions();
mySimulateOptions.rerunExisting = true; 
myWorksheet = simulateWorksheet(myWorksheet, mySimulateOptions);

% Now we create a mapelOptions structure
help mapelOptions
myMapelOptions = mapelOptions;

% First we need to calculate the summary statistics.  In the published
% MAPEL version, we matched mean, standard deviation, and binned data.
myMapelOptions.expData = convertResponseTypeToExpDataTable(myWorksheet, 'N87_agx', 'USUBJID');
myMapelOptions.mnSDTable = convertExpDataToMnSDTable(myMapelOptions);
% Note that this will automatically generate the bin table for now.
% Sometimes, we may want to create our own bin tables rather than using the
% default bin edges.
myMapelOptions.binTable = convertExpDataToBinTable(myMapelOptions);

[nAxis,nVP] = size(myWorksheet.axisProps.axisVP);
% We set additional parameters for the optimization
myMapelOptions.optimizeType = 'pso';
% nVPs is the largest total number of bins we should create.
% We need at least 1 VP per bin
% or we risk assigning weight to areas
% not represented by VPs.  This sets the number of bins on a
% per-axis basis, so there are a total 4 * 13 total population bins
myMapelOptions.nBins = 4;
% We usually don't need very small tolerances to find good solutions, which
% we can assess by the GOF statistic and effN.  The objective is on a
% logarithmic scale.  The optimization may time out before convergence,
% we need to check the GOF to see if the solution is acceptable in terms
% of agreement with data.  Multiple VPop solutions may be acceptable.
myMapelOptions.tol = 1E-2;
myMapelOptions.randomStart = true;
myMapelOptions.spreadOut = 1;
myMapelOptions.minEffN = 30;
myMapelOptions.optimizePopSize = 1000;

% We may want to run a short optimization to create a VPop object
% just to see how well our experimental and simulated data ranges agree
myMapelOptions.optimizeTimeLimit = 10;
myVPop = mapel(myWorksheet, myMapelOptions);

% This diagnostic plot for the population doesn't require a VPop solution, 
% but it does require
% the VPop have both the experimental and simulation data populated.
% ExpData is usually populated manually, but simData is usually populated 
% by the getSimData method of the VPop in the call to mapel.
help plotSimExpDataVPop
plotSimExpDataVPop(myVPop)
% Outputs/interventions/times where the red circles do not cover the
% ranges of the blue triangles well suggest there is not good model/cohort
% coverage and the population fitting have issues.  You need to assess
% whether this means the cohort need to be expanded, whether this is an
% issue with the data (for example, if there is non-mechanistic measurement
% noise in the data that we don't capture), or if there could be an issue
% with the model that we need to fix, discuss with experimenters, or
% accept.  If we decide we want to accept an issue and we do not want it to
% influence fits, we need to set the corresponding weight to 0.


% We had identified several time/experiment datapoints that would be 
% difficult to match due to issues in the experimental data or model 
% behaviors.  We can set the influence of these on the optimization to 0.
myMapelOptions.mnSDTable{1,'weightSD'}=0;
myMapelOptions.mnSDTable{1,'weightMean'}=0;
myMapelOptions.binTable{1,'weight'}=0;

myMapelOptions.mnSDTable{2,'weightSD'}=0;
myMapelOptions.mnSDTable{2,'weightMean'}=0;
myMapelOptions.binTable{2,'weight'}=0;

myMapelOptions.mnSDTable{3,'weightSD'}=0;
myMapelOptions.mnSDTable{3,'weightMean'}=0;
myMapelOptions.binTable{3,'weight'}=0;

myMapelOptions.mnSDTable{26,'weightSD'}=0;
myMapelOptions.mnSDTable{26,'weightMean'}=0;
myMapelOptions.binTable{26,'weight'}=0;

myMapelOptions.mnSDTable{27,'weightSD'}=0;
myMapelOptions.mnSDTable{27,'weightMean'}=0;
myMapelOptions.binTable{27,'weight'}=0;

myMapelOptions.mnSDTable{28,'weightSD'}=0;
myMapelOptions.mnSDTable{28,'weightMean'}=0;
myMapelOptions.binTable{28,'weight'}=0;

% Convergence times may vary between systems
myMapelOptions.optimizeTimeLimit = 2*60*60;
% We can iterate through a few times in case we have a bad seed
% population in the fitting.
for myTestCounter = 1 : 1
        mySaveName = ['example06_',num2str(myTestCounter)];
        myMapelOptions.intSeed = myTestCounter;
        myVPop = mapel(myWorksheet, myMapelOptions); 
        saveVPop(myVPop,mySaveName);
        % We can also display the summary statistic with useEffN
        myVPop.useEffN = true;
        myVPop = evaluateGOF(myVPop);        
        myVPop.mnSDTable{1,'predN'}    
        myVPop.gof
end

% This is an additional diagnostic table for the places where we have
% experimental and simulation data.  Note that these can become
% rather crowded when we have a lot of distinct
% outputs/interventions/times and don't exclude non-weighted times, but
% can be generally useful to get a quick sense of issues with the
% population fitting.
plotMnSDVPop(myVPop)

% There is also a quick plot to visualize how the prevalence weights
% are distributed in the VPop
plotPWHist(myVPop)

% Note that we may find a good VPop and may want to use this as a seed for
% additional VPops, or we may want to pick up with a VPop fitting
% run that was stopped early.
% For example, we can start from a previous VPop as an initial guess,
% and optionally add noise to the initial probabilities
myVPop = loadVPop('example_vpop');
% Set this to 0 if we want to use the exact same starting point
myRandomStart = 0.01;
for myTestCounter = 1 : 1
    myVPop.intSeed = myTestCounter;
    myTestCounter
    newVPop = restartMapel(myVPop, myRandomStart);
    mySaveName = ['example06_rerun_',num2str(myTestCounter)];
    saveVPop(newVPop,mySaveName);
    newVPop.useEffN = true;
    newVPop = evaluateGOF(newVPop);
    newVPop.mnSDTable{1,'predN'}
    newVPop.gof
end