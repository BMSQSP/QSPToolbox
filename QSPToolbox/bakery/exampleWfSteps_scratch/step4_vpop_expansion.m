% Here, we illustrate workflow for creating VPops using MAPEL

% run on blade 30 clusters

clear all;
clc;

% Toolbox initialization
% toolboxPath = 'DRIVE:/DIRs/QSPToolbox';
% addpath(toolboxPath);
initQSPToolbox;

%%
myPath='.../step1/';
myFileName='myWorksheetScreened_AbShedAgxab1petRTlarger';
myWorksheet=loadWorksheetAutoMerge(myFileName, myPath);

updateBounds = {{'axisLesionLinGrowth', [8.3/3	74.7]}, {'axisLesionExpGrowth', [0.0244/3	0.219]},...
    {'axisInitialTumorMass', [0.0763/3	0.24]}};
myWorksheet = resetAxisBounds(myWorksheet, updateBounds);

% This is a sample run with expandVPopEffN
myVPop = loadVPop('example_vpop'); %    %% calibration based on VPop. unweighted the few lines after myMapelOption
% update the VPop, add expDataID to check for the calibration target, for use in linearCalibration and expandVPopEffN
myMapelOptions = initializeVPopPropertiesToOption(myVPop);
myMapelOptions.initialPWs = -1;
myMapelOptions.initialProbs=[];
myMapelOptions.pwStrategy = 'direct';

%% Lu added these two unweight lines on 11/8/2024
%% Lu added these two unweight lines on 11/29/2024
myMapelOptions.mnSDTable{29,'weightSD'}=1;
myMapelOptions.mnSDTable{29,'weightMean'}=1;
myMapelOptions.binTable{29,'weight'}=1;

myMapelOptions.mnSDTable{31,'weightSD'}=1;
myMapelOptions.mnSDTable{31,'weightMean'}=1;
myMapelOptions.binTable{31,'weight'}=1;


myMapelOptions.mnSDTable{9,'weightSD'}=0;
myMapelOptions.mnSDTable{9,'weightMean'}=0;
myMapelOptions.binTable{9,'weight'}=0;

myMapelOptions.mnSDTable{14,'weightSD'}=0;
myMapelOptions.mnSDTable{14,'weightMean'}=0;
myMapelOptions.binTable{14,'weight'}=0;

myMapelOptions.mnSDTable{20,'weightSD'}=0;
myMapelOptions.mnSDTable{20,'weightMean'}=0;
myMapelOptions.binTable{20,'weight'}=0;

myMapelOptions.mnSDTable{25,'weightSD'}=0;
myMapelOptions.mnSDTable{25,'weightMean'}=0;
myMapelOptions.binTable{25,'weight'}=0;

%% LU add stop here

myExpandVPopEffNOptions = expandVPopEffNOptions;
% An identifying suffix
myExpandVPopEffNOptions.suffix = 'burst0p9';
myExpandVPopEffNOptions.wsIterCounter = 0;
% Normally target 150 and often pause at 50
% in initial run.  We'll stop here too in order to
% adjust expansion sized
myExpandVPopEffNOptions.targetEffN = 40; % 40;
%% % added by Lu 20241127
myMapelOptions.minEffN = 0; % myExpandVPopEffNOptions.targetEffN; % added by Lu 20241127
%%
myExpandVPopEffNOptions.maxNewPerIter= 10;
myExpandVPopEffNOptions.expandCohortSize = 5000;
myExpandVPopEffNOptions.effNDelta = 2;
myExpandVPopEffNOptions.minPVal = 0.9; % 0.2; Lu change
myExpandVPopEffNOptions.nTries = 1;
myExpandVPopEffNOptions.nRetries = 1;
myExpandVPopEffNOptions.verbose = true;
myExpandVPopEffNOptions.restartPVal = 1E-4;
myExpandVPopEffNOptions.expandRandomStart = 0;
myExpandVPopEffNOptions.varyMethod = 'gaussian';
myExpandVPopEffNOptions.resampleStd = 0.005;
myExpandVPopEffNOptions.expandEdgeVPs=true;
myExpandVPopEffNOptions.maxNewPerOld = 3;
% Here, we enforce any checks on sampled VPs,
% before simulation and adjust VPs if needed.
myExpandVPopEffNOptions.screenFunctionName = '';
myExpandVPopEffNOptions.linearExpandFlag = true;
myExpandVPopEffNOptions.minEffNlinearflag = true;
myExpandVPopEffNOptions.linearExpandCohortSize = 5000;
myExpandVPopEffNOptions.maxIterlinearExpand = 50;
 myExpandVPopEffNOptions.nVPMax = 1000; % 1050 in earlier try
 myExpandVPopEffNOptions.minPVallinear = 1e-2;
 myExpandVPopEffNOptions.nCluster = 400;
% This takes a while to run.  Every iteration will
% write to file when a better VPop is found
% and a new worksheet will write to file
% when new VPs are added.
% You can delete intermediate iterations between your start file and
% stopping point.  They are included to allow you to more easily
% track progress.

hotStartVPop='';
[newWorksheet, newVPop] = expandVPopEffN(myWorksheet,myExpandVPopEffNOptions,myMapelOptions,hotStartVPop);

