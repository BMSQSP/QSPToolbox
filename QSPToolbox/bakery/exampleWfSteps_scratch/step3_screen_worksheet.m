% loadWorksheet(myFileName, myPath, myFormat, opts)

clear all;
clc;

myPath='/.../step1/';

% Here we give an example for screening a worksheet
% against responsetypes
% The script is broken into several steps:
%  #1: screening the bulk simulated VPs from multiple worksheets
%       (files not included)
%  #2: clustering
%  #3: finding missing phenotypes (several iterations here)
%       and re-clustering
% Step 3 is somewhat optional.  Sometimes we get good calibrations
% without it.

% Toolbox initialization
initQSPToolbox;

% Info to read in worksheets.  Here we assume the 30 worksheets
% from step 1 have been simulated and have results attached.
% This code is dependent on the previous step, or could be modified for
% other projects.
myFilePrefix = 'myWorksheet';
nWorksheets = 40;

% Try reading in the first one
myWorksheet = loadWorksheet([myFilePrefix,'_1']);
% Alternately, you can read in the sample worksheet just
% for a worksheet with responseTypes attached
% myWorksheet = loadWorksheet(myFileNameSample);

% Set up the bounds on the responseTypes for screening
% the split-up submissions
% Initially we only keep VPs in range.
allResponseTypeIDs = getResponseTypeIDs(myWorksheet);
myScreenTable = createScreenTable(myWorksheet, allResponseTypeIDs, false);

% Now read in the worksheets
% and screen the VPs in a PARFOR
% myWorksheetCell = cell(1,nWorksheets);
for i=1:nWorksheets
    i
    myFileName=['myWorksheet_' num2str(i)];
    myWorksheet=loadWorksheet(myFileName, myPath);    
        myResponseSummaryTable = createResponseSummaryTable(myWorksheet, 'N87_agx');
        myIndex1 = find(ismember(myResponseSummaryTable.rowNames,'culture_internalized'));
        myIndices1 = find(myResponseSummaryTable.values(myIndex1,:)<=1.0);
%          myIndex2 = find(ismember(myResponseSummaryTable.rowNames,'buffer_injection_1200147_growth_assay'));
%          myIndices2 = find(myResponseSummaryTable.values(myIndex2,:)<=1.0);
        myIndex3 = find(ismember(myResponseSummaryTable.rowNames,'buffer_injection_1200147_shed_assay'));
        myIndices3 = find(myResponseSummaryTable.values(myIndex3,:)<=2.0);  %0.75);
%         myIndex4 = find(ismember(myResponseSummaryTable.rowNames,'antibody_injection_1200147_growth_assay'));
%         myIndices4 = find(myResponseSummaryTable.values(myIndex4,:)<=1.0);
%         myIndex5 = find(ismember(myResponseSummaryTable.rowNames,'antibody_injection_1200147_shed_assay'));
%         myIndices5 = find(myResponseSummaryTable.values(myIndex5,:)<=0.75);  %% myIndices5=empty % remove after removing this response type
        myIndex6 = find(ismember(myResponseSummaryTable.rowNames,'agxab1pet_injection'));
        myIndices6 = find(myResponseSummaryTable.values(myIndex6,:)<=3.0);
    %     myIndices = intersect(intersect(intersect(intersect(intersect(myIndices1,myIndices2),myIndices3),myIndices4),myIndices5),myIndices6);
%         myIndices = intersect(intersect(intersect(myIndices1,myIndices3),myIndices5),myIndices6);
%           myIndices = intersect(intersect(myIndices1,myIndices3),myIndices6);
         myIndices = intersect(intersect(myIndices1,myIndices3),myIndices6);

        length(myIndices);
        allVPIDs = getVPIDs(myWorksheet);
        myVPIDs = allVPIDs(myIndices);
        myWorksheet = copyWorksheet(myWorksheet, myVPIDs);
        length(myVPIDs)
        if i==1
            myWorksheet_all = myWorksheet;
        else
            myWorksheet_all = mergeWorksheets(myWorksheet_all,myWorksheet);
        end
end

saveWorksheetAutoSplit(myWorksheet_all,'myWorksheetScreened_AbShedAgxab1petRTlarger');









