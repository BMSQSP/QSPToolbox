% modified the example04 script ADC model to step1 based on the reference workflow
% use example06 worksheet's setup
clear all;
clc;
% Toolbox initialization
% toolboxPath = 'DRIVE:/DIRs/QSPToolbox';
% addpath(toolboxPath);
initQSPToolbox;

% Read in the model
myWorksheet = readQSPModel('adc_platform_public_v0p15p1.sbproj');
% Create initial blank VP
% myWorksheet = createVPs(myWorksheet,{'initialVP'},{{}}); % if use this will end up NaN 
myWorksheet = readVPTable(myWorksheet, 'example_vp_file.txt'); % change to example02 lines to try
% Read in the intervention (row) definitions
myWorksheet = readInterventionTable(myWorksheet, 'example_intervention_file.txt');

% Set some of the simulation properties - how often to save a result,
% how long to simulate for

% add in target time point to ensure we store the simulation results at the correct time points
myMapelOptions = mapelOptions;
myWorksheetFileName = 'example_cohort_worksheet1'; % check this RT?
exampleWorksheet = loadWorksheet(myWorksheetFileName);
myMapelOptions.expData = convertResponseTypeToExpDataTable(exampleWorksheet, 'N87_agx', 'USUBJID');
myWorksheet.simProps.sampleTimes = sort(unique([0:.025:35,myMapelOptions.expData.time'])','ascend');

% We can also adjust tolerances as needed
mySimulateOptions = simulateOptions();
mySimulateOptions.rerunExisting = true; 
myWorksheet.simProps.relativeTolerance = 1E-6;

%% Set mech. axes
% Add Axis Definition
AxisRange = readtable('axes.txt','TreatAsEmpty',{'.'});
for i = 1:size(AxisRange,1)
    myWorksheet = addAxisDef(myWorksheet, AxisRange.axisID{i}, ...
        {AxisRange.elementNames{i}, AxisRange.elementTypes{i}, ...
        [AxisRange.bounds1_1(i), AxisRange.bounds1_2(i)], AxisRange.scale1{i}});
end

% % Add mechanistic axes definitions to the worksheet.
% % Note that adding axes will trigger a model export and getting the model
% % "elements" (parameters, species, compartment values)
% myWorksheet = addAxisDef(myWorksheet,'axisShedding',{'par_ag_shed_rate_t','parameter', [0, 10]});
% myWorksheet = addAxisDef(myWorksheet,'axisInternalization',{'par_ag_endocytosis','parameter', [0, 10]});
% myWorksheet = addAxisDef(myWorksheet,'axisWallPoreB',{'par_fraction_capillary_pore_area_per_unit_thickness_b','parameter', [0.2, 0.8]});
% myWorksheet = addAxisDef(myWorksheet,'axisShedClearance',{'par_shed_antigen_clearance_p1','parameter', [20, 80]});
% myWorksheet = addAxisDef(myWorksheet,'axisLabelExport',{'par_payload_secretion_rate','parameter', [0, 24]});

%% Read experiment data
% Add experimental datasets to the worksheet, so we can include
% them on the plot and include them in response types.
myWorksheet = readExperimentData(myWorksheet, 'fluorescence_quench_internalization_data_n87_ab1.txt');
myWorksheet = readExperimentData(myWorksheet, 'xenograft_data_1200115_n87_buffer_only.txt');
myWorksheet = readExperimentData(myWorksheet, 'xenograft_data_1200147_n87_buffer_only.txt');
myWorksheet = readExperimentData(myWorksheet, 'xenograft_data_1200147_n87_ab1_only.txt');
myWorksheet = readExperimentData(myWorksheet, 'xenograft_data_95305073_n87.txt');
% You can see the experimental data in the worksheet, it is read as tables.
% This does create some name length restrictions for the experimental
% variables
myWorksheet.expData{:}

%% Set response types
% Response types serve as a mechanism to score how well we agree with data.
% We add one here.
myWorksheet.responseTypes={};
myWorksheet = createResponseType(myWorksheet,'N87_agx');
PointsRTRange = readtable('pointsRT.txt','TreatAsEmpty',{'.'});
for i=1:size(PointsRTRange,1)
   RTi=PointsRTRange(i,:)
   myWorksheet = addResponseTypeElement('points',myWorksheet,...
       {RTi.rtID{1},'N87_agx',RTi.modelYVar{1},RTi.modelYVarType{1},...
       RTi.interventionID{1},RTi.expDataID{1},RTi.expDataTimeVar{1},...
       RTi.expDataYVar{1},RTi.objectiveType{1},RTi.weight}); 
end

% Lu: to remove bound responseTypes, for example: antibody_injection_1200147_shed_assay because 2 time/data points always are out of bounds
idx=find(ismember(getResponseTypeElementIDs(myWorksheet,'N87_agx'),'antibody_injection_1200147_shed_assay'));
myWorksheet.responseTypes{1}.elements(idx) = []; 
% saveWorksheetAutoSplit(myWorksheet,'example_cohort_worksheet1_rvAbInjection_result');

% Check the response type to make sure we have paired
% everything correctly
myPassCheck = verifyResponseType(myWorksheet, 'N87_agx');

%% burst simulation
saveIDlist = {};
myVaryAxesOptions = varyAxesOptions;
myVaryAxesOptions.baseVPIDs = getVPIDs(myWorksheet);
myVaryAxesOptions.newPerOld = 2e5;

myVaryAxesOptions.intSeed = 1;
myVaryAxesOptions.varyMethod = 'sobol';
myVaryAxesOptions.varyAxisIDs = getAxisDefIDs(myWorksheet);
myVaryAxesOptions.baseVPIDs = {'n87_xenograft_scenario_0006'}; % added as example02. try?
myWorksheet = addVariedVPs(myWorksheet,myVaryAxesOptions);
allVPIDs = getVPIDs(myWorksheet);
myWorksheet = copyWorksheet(myWorksheet,allVPIDs(2:end));
% Confirm between 0 and 1
min(myWorksheet.axisProps.axisVP.coefficients')
max(myWorksheet.axisProps.axisVP.coefficients')
% Make sure we're not missing variables that are clearly needed
% and that there are no duplicates
allResponseTypeElements = getResponseTypeElementIDs(myWorksheet,'N87_agx');
for elementCounter = 1:length(allResponseTypeElements)
    saveIDlist = [saveIDlist,myWorksheet.responseTypes{1}.elements{elementCounter}.modelYVar];
end
myWorksheet.simProps.saveElementResultIDs=unique([myWorksheet.simProps.saveElementResultIDs,saveIDlist],'stable');
% Accelerate the model before saving so speed simulation steps.
% Note this is processor/platform specific, so if you run on linux
% to save time you should compile/accelerate on linux.
help compileModel
myWorksheet = compileModel(myWorksheet, true)
saveWorksheet(myWorksheet,'myWorksheet_prescreened');

%%
% Now get ready to run, here we split into batches and
% you can run alternate worksheets on different servers.
% Hopefully you have a shared drive to read from
initQSPToolbox
% Break the worksheet into smaller sets of 10,000 VPs for simulations.
% This just helps in case there are failed runs and to reduce input-output times.
myWorksheet = loadWorksheet('myWorksheet_prescreened');
%simRunVPSize = 10000;
simRunVPSize = 5000;
allVPIDs = getVPIDs(myWorksheet);
myWorksheet=compileModel(myWorksheet);
nVPs = length(getVPIDs(myWorksheet));
nWorksheets = ceil(nVPs/simRunVPSize)
for splitCounter = 1 : nWorksheets
    curVPIDs = allVPIDs(((splitCounter-1)*simRunVPSize+1):((splitCounter)*simRunVPSize));
    simWorksheet = copyWorksheet(myWorksheet, curVPIDs, false, false);
    saveWorksheet(simWorksheet,['myWorksheet_',num2str(splitCounter)]);
end

%%
% Simulate in batches.
% Each worksheet takes about 30 minutes to run on a 64-core
% server with Xeon E5-6660 v4 processors with v3b4mdsc of the model
% The next steps will involve loading the worksheets, screening VPs,
% adding any missing phenotypes, and then VPop expansions
% SERVER 1
initQSPToolbox
mySimulateOptions=simulateOptions;
mySimulateOptions.poolClose = false;
mySimulateOptions.poolRestart = false;
for splitCounter = 1:nWorksheets
    simWorksheet = loadWorksheet(['myWorksheet_',num2str(splitCounter)]);
    tic
    simWorksheet = simulateWorksheet(simWorksheet,mySimulateOptions);
    toc
    saveWorksheet(simWorksheet,['myWorksheet_',num2str(splitCounter)]);
end
