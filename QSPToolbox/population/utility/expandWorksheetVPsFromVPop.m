function [myWorksheet, newPassNames] = expandWorksheetVPsFromVPop(myWorksheet,newVPop, myMapelOptions,suffix,wsIterCounter, maxNewPerIter, testBounds, expandCohortSize, varyMethod, gaussianStd, maxNewPerOld, nUnweightedParents, selectByParent, myScreenFunctionName)
% This function expands a worksheet given a VPop.  It selected out VPs to expand around,
% samples for new VPs, scores the "children" based on available data, and adds
% the best children to the worksheet.
%
%  myWorksheet
%  newVPop
%  myMapelOptions
%  suffix:             a text descriptor string that will be included in what 
%                       will be written to file.  This is also used
%                       in setting VP identities.
%  wsIterCounter:      tracks the iterations through the algorithm.  Keep
%                       incrementing to avoid issues with repeated VPIDs.
%  maxNewPerIter:      maximum new VPs we can add per iteration.  Set to 
%                       -1 to use the VPop effN
%  testBounds:         a cell array of bounds for the response Types.
%                       there's a quirky formatting here where each cell
%                       is the limits on the coutput of non-axis
%                       rows of the values field of createResponseSummaryTable.
%                       TODO: update this to take a cell array of
%                       standard outputs of evaluateResponseType
%  expandCohortSize:   size of the cohort to generate for testing
%  varyMethod:         method for resampling.  i.e. 'gaussian' or 'localPCA'
%  gaussianStd:        standard deviation for the re-sampled parameters.
%                       note this is applied across all axes in the transformed
%                       units (i.e. within bounds normalized 0-1).
%  maxNewPerOld:       maximum number of children per weighted parent
%  nUnweightedParents: Number of unweighted VPs to try as a parent 
%                       for expansion each iteration.
%  selectByParent:     a boolean variable indicating whether children
%                       VPs are selected for inclusion each iteraction 
%                       based on their parent, or just pooled together
%  myScreenFunctionName: a string indicating a function to use for screening
%                         VPs before simulation.  It should take two input
%                         arguments: a worksheet and a list of VPIDs to
%                         screen.  It should return a worksheet with
%                         identical number of names of VPs, possibly 
%                         modified after bing checked against some
%                         criteria. '' indicates no screening.
%                      
%                      
% RETURNS:
%  myWorksheet
%  newPassNames
%
% TODO: this function could use more input proofing or an options object.  This
%       is also mainly intended to be called from other functions, so it hasn't
%       been done yet.

originalVPIDs = getVPIDs(myWorksheet); 
if length(originalVPIDs) > length(newVPop.pws)
    warning(['More VPs in worksheet than PWs in VPop in ',mfilename,'.  Proceeding assuming the PWs correspond to the first VPs in the worksheet.'])
    continueFlag = true;
elseif length(originalVPIDs) < length(newVPop.pws)
    warning(['Fewer VPs in worksheet than PWs in VPop in ',mfilename,'.  Exiting...'])
    continueFlag = false; 
else
    continueFlag = true;
end

    
if continueFlag    
    
    myCoeffs = getVPCoeffs(myWorksheet);
    [nAxis, nOriginalVPs] = size(myCoeffs);	
    
    % We will pre-rank parent VPs we might want to add
    % to resampling based on the observed data    
    % This step will only be impactful if
    % nUnweightedParents > 0
    % We need to create an updated the VPop in case
    % new VPs were added to the worksheet
    if isa(newVPop,'VPop')
        testVPop = VPop;
    elseif isa(newVPop,'VPopRECIST')
        testVPop = VPopRECIST;
    elseif isa(newVPop,'VPopRECISTnoBin')
        testVPop = VPopRECISTnoBin;
    end
    testVPop.expData = newVPop.expData;
    testVPop.mnSDTable = newVPop.mnSDTable;
    testVPop.binTable = newVPop.binTable;
    testVPop.distTable = newVPop.distTable;
    testVPop.distTable2D = newVPop.distTable2D;   
    testVPop.corTable = newVPop.corTable;
    testVPop.subpopTable = newVPop.subpopTable;
    if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
        testVPop.brTableRECIST = newVPop.brTableRECIST;
        testVPop.rTableRECIST = newVPop.rTableRECIST;        
        testVPop.relSLDvar = newVPop.relSLDvar;
        testVPop.absALDVar = newVPop.absALDVar;
        testVPop.crCutoff = newVPop.crCutoff;             
        testVPop.recistSimFilter = createRECISTSimFilter(myWorksheet, testVPop);
    end
    if ~isa(newVPop,'VPopRECISTnoBin')
		testVPop = testVPop.assignIndices(myWorksheet, myMapelOptions);
    end
    testVPop = testVPop.getSimData(myWorksheet);
	testVPop.subpopTable = updateSubpopTableVPs(testVPop.subpopTable, myWorksheet);
    testVPop = testVPop.addTableSimVals();  
    % For evaluation
    % coerce pws to be the same.
    testVPop.pws = ones(1,nOriginalVPs)./nOriginalVPs;
    testVPop = testVPop.addPredTableVals();    
    
    % Get the scores
    originalVPScores = scoreWorksheetVPs(testVPop,1:nOriginalVPs,1:nOriginalVPs);
    [nScoresPerVP, ~] = size(originalVPScores);
    
    % Rank the original VPs by score
    sortIndices = nan(nScoresPerVP, nOriginalVPs);
    for rowCounter = 1 : nScoresPerVP
        nonZeroIndices = find(originalVPScores(rowCounter,:)>0);
        [curScores, curIndices] = sort(originalVPScores(rowCounter,nonZeroIndices),'descend');
        originalVPScores(rowCounter,:) = nan(1,nOriginalVPs);
        originalVPScores(rowCounter,1:length(nonZeroIndices)) = curScores;
        sortIndices(rowCounter,1:length(nonZeroIndices)) = nonZeroIndices(curIndices);
    end
    sortIndices=reshape(sortIndices,1,[]);
    sortIndices = sortIndices(find(~isnan(sortIndices)));
    sortIndices = unique(sortIndices,'stable');
    originalVPIDsSort = originalVPIDs(sortIndices);
    
    % Consider VPs for inclusion as seed if they are weighted "heavily"
    pwExpandCutoff = 0.01;
    allResponseTypeIDs = getResponseTypeIDs(myWorksheet);
    nResponseTypes = length(allResponseTypeIDs);
    myPWs = newVPop.pws;
    curEffN = round(1/sum(myPWs.^2));
    highVPindices1 = find(myPWs>=pwExpandCutoff);
    [~, highVPindices2] = sort(myPWs, 'descend');
    highVPindices2 = highVPindices2(1 : curEffN);
    highVPindicesCat = [highVPindices2,highVPindices1];
    [highVPindices,i,j] = unique(highVPindicesCat, 'first'); 
    highVPindices = highVPindicesCat(sort(i));
    highVPIDs = originalVPIDs(highVPindices);
    [originalVPIDsSort,sortIndicesPick] = setdiff(originalVPIDsSort,highVPIDs,'stable');
    sortIndices = sortIndices(sortIndicesPick);
    % Now combine considerations for heavily weighted VPs and
    % parents that look like they are useful.
    highVPIDs = [highVPIDs,originalVPIDsSort(1:min(nUnweightedParents,length(originalVPIDsSort)))];
    highVPindices = [highVPindices,sortIndices(1:min(nUnweightedParents,length(originalVPIDsSort)))];

    % Decide how many of the VPs from the seed expansion we can add.
    if maxNewPerIter < 0
        maxNewPerIterChecked = curEffN;
    else
        maxNewPerIterChecked = maxNewPerIter;
    end
	
    % We'll allow more than maxNewPerOld if it looks like we want
    % to allow many new VPs per seed.  We won't allow more
	% if we are selecting children based on parent and
	% we have many parents. If we are pooling all the children
	% together and picking the best, we will ignore this.
    maxNewPerOld = max(maxNewPerOld,ceil(maxNewPerIterChecked/length(highVPindices)));

    myVaryAxesOptions = varyAxesOptions;
    myVaryAxesOptions.varyMethod = varyMethod;
    myVaryAxesOptions.gaussianStd = gaussianStd;
    myVaryAxesOptions.varyAxisIDs = getAxisDefIDs(myWorksheet);
    myVaryAxesOptions.intSeed = wsIterCounter;

    myVaryAxesOptions.additionalIDString = [suffix,num2str(wsIterCounter)];
    myVaryAxesOptions.baseVPIDs = highVPIDs;
    myVaryAxesOptions.newPerOld = ceil(expandCohortSize/length(highVPIDs));

    jitteredWorksheet = addVariedVPs(myWorksheet, myVaryAxesOptions);
    curVPIDs = getVPIDs(jitteredWorksheet);
    newIndices = (find(~ismember(curVPIDs,originalVPIDs)));
    newVPIDs = curVPIDs(newIndices);

    % We will also randomize a coefficient
    for vpCounter = 1 : length(newVPIDs);
        curVPID = newVPIDs{vpCounter};
        % randomize one of the new VP axis coefficients
        axisIndex = randsample([1:nAxis],1);
        vpIndex = find(ismember(curVPIDs,curVPID));
        jitteredWorksheet.axisProps.axisVP.coefficients(axisIndex,vpIndex) = rand(1);
    end

    mySimulateOptions = simulateOptions;
    mySimulateOptions.rerunExisting = false;
    mySimulateOptions.optimizeType = 'none';
	% Inherit the pool properties
	mySimulateOptions.poolRestart = myMapelOptions.poolRestart;
	mySimulateOptions.poolClose = myMapelOptions.poolClose;    
    
    % Also screen the worksheet if a function is provided
    if length(myScreenFunctionName) > 0
        jitteredWorksheet = eval([myScreenFunctionName,'(jitteredWorksheet,newVPIDs,mySimulateOptions)']);
    end
    
    jitteredWorksheet = simulateWorksheet(jitteredWorksheet,mySimulateOptions);

    % Update VPIDs as some will fail simulation
    curVPIDs = getVPIDs(jitteredWorksheet);
    newIndices = (find(~ismember(curVPIDs,originalVPIDs)));
    newVPIDs = curVPIDs(newIndices);

    % Identify the VPs that don't fulfill the worksheet response
    % as well as those in the initial worksheet in order to filter
    % them.
    newInvalidIndices = nan(1,0);
    for responseTypeCounter = 1 : nResponseTypes
        curTable = createResponseSummaryTable(jitteredWorksheet, allResponseTypeIDs{responseTypeCounter});
        curInvalidIndices = find(sum(curTable.values((nAxis+1):end,newIndices)>testBounds{responseTypeCounter},1)>0);
        if length(curInvalidIndices) > 0
            newInvalidIndices = [newInvalidIndices, curInvalidIndices];
        end
    end

    % Get the new VPs in the worksheet that look OK
    allChildrenNames = newVPIDs;
    if length(newInvalidIndices) > 0
        newInvalidIndices = sort(unique(newInvalidIndices),'ascend');
        newInvalidIDs = newVPIDs(newInvalidIndices);
        newValidIndices = find(~ismember(allChildrenNames,newInvalidIDs));
        allChildrenNames = allChildrenNames(newValidIndices);
    else
        newValidIndices = nan(1,0);
        allChildrenNames = cell(1,0);
    end

    % Just keep the valid VPs
    jitteredWorksheet = copyWorksheet(jitteredWorksheet,[originalVPIDs,allChildrenNames]);
    curVPIDs = getVPIDs(jitteredWorksheet);
    originalIndices = (find(ismember(curVPIDs,originalVPIDs)));
    newIndices = (find(~ismember(curVPIDs,originalVPIDs)));
    newVPIDs = curVPIDs(newIndices);

    
    
    
    % We create a "dummy" vpop object to help identify
    % areas where we may need to expand the range in simulated
    % outcomes
    if isa(newVPop,'VPop')
        testVPop = VPop;
    elseif isa(newVPop,'VPopRECIST')
        testVPop = VPopRECIST;
    elseif isa(newVPop,'VPopRECISTnoBin')
        testVPop = VPopRECISTnoBin;
    end
    testVPop.expData = newVPop.expData;
    testVPop.mnSDTable = newVPop.mnSDTable;
    testVPop.binTable = newVPop.binTable;
    testVPop.distTable = newVPop.distTable;
    testVPop.distTable2D = newVPop.distTable2D;   
    testVPop.corTable = newVPop.corTable;
    testVPop.subpopTable = newVPop.subpopTable;    
    if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
        testVPop.brTableRECIST = newVPop.brTableRECIST;
        testVPop.rTableRECIST = newVPop.rTableRECIST;        
        testVPop.relSLDvar = newVPop.relSLDvar;
        testVPop.absALDVar = newVPop.absALDVar;
        testVPop.crCutoff = newVPop.crCutoff;             
        testVPop.recistSimFilter = createRECISTSimFilter(jitteredWorksheet, testVPop);
    end
    if ~isa(newVPop,'VPopRECISTnoBin')
		testVPop = testVPop.assignIndices(jitteredWorksheet, myMapelOptions);
    end
    testVPop = testVPop.getSimData(jitteredWorksheet);
	testVPop.subpopTable = updateSubpopTableVPs(testVPop.subpopTable, jitteredWorksheet);
    testVPop = testVPop.addTableSimVals();  
    % For evaluation
    % coerce pws to be the same.
    testVPop.pws = ones(1,length(curVPIDs))./length(curVPIDs);
    testVPop = testVPop.addPredTableVals();    
    
	% Now get the score matrix for the new VPs
    % We can remove the LC scoring option for now
    % if isa(newVPop,'VPopRECIST')
    %     newVPScores = scoreWorksheetVPsLC(testVPop,originalIndices,newIndices);
    % else
        newVPScores = scoreWorksheetVPs(testVPop,originalIndices,newIndices);
    % end


    newPassNames = cell(1,0);
	if selectByParent
		% Select valid VPs from each higher weighted seed VP
		for highCounter = 1 : length(highVPIDs)
			parentID=highVPIDs(highCounter);
			childrenBase = strcat(parentID, ['_',suffix,num2str(wsIterCounter)]);
			allChildrenBaseIndices = cellfun(@isempty,strfind(curVPIDs,childrenBase));
			allChildrenBaseIndices=find(~allChildrenBaseIndices);
			allChildrenBaseIndices=intersect(allChildrenBaseIndices,newIndices);
			childIDs = curVPIDs(allChildrenBaseIndices);
			childIDs = newVPIDs(find(ismember(newVPIDs,childIDs)));
			childScores = newVPScores(:,find(ismember(newVPIDs,childIDs)));
			% We will prioritize VPs that score well
			[nScoresPerVP, nCurChildren] = size(childScores);
			sortIndices = nan(nScoresPerVP, nCurChildren);
			for rowCounter = 1 : nScoresPerVP
				nonZeroIndices = find(childScores(rowCounter,:)>0);
				[curScores, curIndices] = sort(childScores(rowCounter,nonZeroIndices),'descend');
				childScores(rowCounter,:) = nan(1,nCurChildren);
				childScores(rowCounter,1:length(nonZeroIndices)) = curScores;
				sortIndices(rowCounter,1:length(nonZeroIndices)) = nonZeroIndices(curIndices);
			end
			sortIndices=reshape(sortIndices,1,[]);
			sortIndices = sortIndices(find(~isnan(sortIndices)));
			sortIndices = unique(sortIndices,'stable');        
			childIDs = childIDs(sortIndices);
			npass = length(childIDs);
			if (npass > 0)
				newPassNames = [newPassNames,childIDs(1:min(npass,maxNewPerOld))];
			end
		end
		if maxNewPerIterChecked < 1
			newPassNames = cell(1,0);
		elseif length(newPassNames) > maxNewPerIterChecked
			newPassNames = newPassNames(1 : maxNewPerIterChecked);
		elseif length(newPassNames) < maxNewPerIterChecked
			addIndices = find(~ismember(newVPIDs,newPassNames));
			if length(addIndices) > 0
				addVPIDs = newVPIDs(addIndices);
				addScores = newVPScores(addIndices);
				[addScores, addIndices] = sort(addScores, 'descend');
				nToAdd = min(maxNewPerIterChecked-length(newPassNames),length(addScores));
				addVPIDs = addVPIDs(addIndices(1:nToAdd));
				newPassNames = [newPassNames, addVPIDs];
			end
		end		
	else
		% Otherwise, we pool the children from
		% all of our test VPs and just take the best
		[nScoresPerVP, nCurChildren] = size(newVPScores);
		childScores = newVPScores;
		sortIndices = nan(nScoresPerVP, nCurChildren);
		for rowCounter = 1 : nScoresPerVP
			nonZeroIndices = find(childScores(rowCounter,:)>0);
			[curScores, curIndices] = sort(childScores(rowCounter,nonZeroIndices),'descend');
			childScores(rowCounter,:) = nan(1,nCurChildren);
			childScores(rowCounter,1:length(nonZeroIndices)) = curScores;
			sortIndices(rowCounter,1:length(nonZeroIndices)) = nonZeroIndices(curIndices);
		end
		sortIndices=reshape(sortIndices,1,[]);
		sortIndices = sortIndices(find(~isnan(sortIndices)));
		sortIndices = unique(sortIndices,'stable');        
		childIDs = newVPIDs(sortIndices);
		npass = length(childIDs);
		% We just allow up to maxNewPerIter
		% in the case where we are pooling children
		if (npass > 0)
			newPassNames = childIDs(1:min(npass,maxNewPerIter));
		end
	end

    % The repeats should already be screened, but just in case
    newPassNames = unique(newPassNames,'first');

    % Finalize the updated cohort
    mergeNames = [originalVPIDs,newPassNames];
    myWorksheet = copyWorksheet(jitteredWorksheet,mergeNames);

    %         % Repopulate results.  Note that this is not necessary
    %         % and some time to each iteration, but help to avoid issues
    %         % if the original worksheet was carrying old, invalid results.
    %         % This is therefore removed, better to check the worksheet
    %         % VPs before calling expandVPopEffN.
    %         myWorksheet.results = {};
    %         myWorksheet = simulateWorksheet(myWorksheet);
else
    warning(['Unable to proceed in ',mfilename,'.  Returning input worksheet.'])
end
end