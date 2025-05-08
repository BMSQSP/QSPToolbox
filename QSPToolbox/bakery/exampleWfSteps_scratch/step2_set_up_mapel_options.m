% Toolbox initialization
initQSPToolbox;

addpath("../step1")


myWorksheet = loadWorksheet('myWorksheet_1');

myMapelOptions = mapelOptions;

% First we need to calculate the summary statistics.  In the published
% MAPEL version, we matched mean, standard deviation, and binned data.
myMapelOptions.expData = convertResponseTypeToExpDataTable(myWorksheet, 'N87_agx', 'USUBJID');
myMapelOptions.mnSDTable = convertExpDataToMnSDTable(myMapelOptions);
% We create a default subpopTable for 'all', since there are noise
% subpopulations of special interest for this calibration.
myMapelOptions.subpopTable = createSubpopTable(myWorksheet);
% Note that this will automatically generate the bin table for now.
% Sometimes, we may want to create our own bin tables rather than using the
% default bin edges.
% Note the bin functions have changes since Cheng et al. AAPS J 2017,
% so bin edges may vary slightly from the available historical
% VPOP used later in this example.
myMapelOptions.pwStrategy = 'direct'; %
myMapelOptions.binTable = convertExpDataToBinTable(myMapelOptions);
[nAxis,nVP] = size(myWorksheet.axisProps.axisVP);
% nVPs is the largest total number of bins we should create.
% We need at least 1 VP per bin
% or we risk assigning weight to areas
% not represented by VPs.  This sets the number of bins on a
% per-axis basis, so there are a total 4 * 13 total population bins
myMapelOptions.nBins = 4;
% We set additional parameters for the optimization
myMapelOptions.optimizeType = 'pso';
% We usually don't need very small tolerances to find good solutions, which
% we can assess by the GOF statistic and effN.  The objective is on a
% logarithmic scale.  The optimization may time out before convergence,
% we need to check the GOF to see if the solution is acceptable in terms
% of agreement with data.  Multiple VPop solutions may be acceptable.
myMapelOptions.tol = 1E-2;
myMapelOptions.randomStart = true;
myMapelOptions.spreadOut = 1;
myMapelOptions.minEffN = 0; % 30
myMapelOptions.optimizePopSize = 1000;

% We may want to run a short optimization to create a VPop object
% just to see how well our experimental and simulated data ranges agree
myMapelOptions.optimizeTimeLimit = 10;

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

% Save the results to file:
saveMapelOptions(myMapelOptions, 'myMapelOptions')

