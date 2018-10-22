function [myWorksheet, newPassNames] = expandWorksheetVPsFromVPopTest(myWorksheet,newVPop, myMapelOptions,suffix,wsIterCounter, maxNewPerIter, testBounds, expandCohortSize, gaussianStd, maxNewPerOld)
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
%  gaussianStd:        standard deviation for the guassian sampled parameters.
%                       note this is applied across all axes in the transformed
%                       units (i.e. within bounds normalized 0-1).
%  maxNewPerOld:       maximum number of children per weighted parent
%                      
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
    % Consider VPs for inclusion as seed if they are weighted this heavily
    pwExpandCutoff = 0.01;

    allResponseTypeIDs = getResponseTypeIDs(myWorksheet);
    nResponseTypes = length(allResponseTypeIDs);
    myCoeffs = getVPCoeffs(myWorksheet);
    [nAxis, ~] = size(myCoeffs);
    myPWs = newVPop.pws;
    curEffN = round(1/sum(myPWs.^2));
    highVPindices1 = find(myPWs>=pwExpandCutoff);
    [~, highVPindices2] = sort(myPWs, 'descend');
    highVPindices2 = highVPindices2(1 : curEffN);
    highVPindicesCat = [highVPindices2,highVPindices1];
    [highVPindices,i,j] = unique(highVPindicesCat, 'first'); 
    highVPindices = highVPindicesCat(sort(i));
    highVPIDs = originalVPIDs(highVPindices);

    % Decide how many of the VPs from the seed expansion we can add.
    if maxNewPerIter < 0
        maxNewPerIterChecked = curEffN;
    else
        maxNewPerIterChecked = maxNewPerIter;
    end

    % Maximum target number of VPs resulting from each seed to keep.
	% For RECIST, we also check VPs mechanistically dropping out (RSCORE)
	% Several strategies were tried here:
	% V792 we have higher newPerOld
    % if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
        % maxNewPerOld = 4;
    % else
		% maxNewPerOld = 3;
	% end
	% % For V793 try limiting the min maxNewPerOld to just 2
	% maxNewPerOld = 2;
	% For V794 try lowering the min maxNewPerOld even more to just 1
	% maxNewPerOld = 1;	
	% For V795 we go back to to strategy of V793 (2 children per parent max)
	% maxNewPerOld = 2;
	% V801 we have higher newPerOld
    % if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
    %   maxNewPerOld = 3;
    % else
    %   maxNewPerOld = 2;
    % end    
    % For V802 we revert (2 children per parent max).  This actually looks
    % like it was working well.
	% maxNewPerOld = 2;
    % For V811 we revert (3 children per parent max).  This 
    % is to allow for inclusion of 2D evaluations.
	% For V831 we take this as an input argument.
	% maxNewPerOld = 3	
	
    % We'll allow more if it looks like we want to allow many new VPs per seed.
    maxNewPerOld = max(maxNewPerOld,ceil(maxNewPerIterChecked/length(highVPindices)));

    myVaryAxesOptions = varyAxesOptions;
    myVaryAxesOptions.varyMethod = 'gaussian';
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
    jitteredWorksheet = simulateWorksheet(jitteredWorksheet,mySimulateOptions);

    % Update VPIDs as some will fail simulation
    curVPIDs = getVPIDs(jitteredWorksheet);
    newIndices = (find(~ismember(curVPIDs,originalVPIDs)));
    newVPIDs = curVPIDs(newIndices);

    % identify the VPs that don't fulfill the worksheet response
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
    newInvalidIndices = nan(1,0);

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
    if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
        testVPop.distTable2D = newVPop.distTable2D;
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
    testVPop.pws = ones(1,length(curVPIDs))./length(curVPIDs);
    testVPop = testVPop.addDistTableSimVals();  
	flagCheck2D = false;
	if ~isa(testVPop,'VPop')
        testVPop = testVPop.addDistTable2DSimVals();
		flagCheck2D = true;
    end
    testVPop = testVPop.addPredTableVals();    

    [nTestOutcomes,~] = size(newVPop.simData.rowInfo);
    ExpandIncreaseSimDataRows = nan(0,1);
    ExpandDecreaseSimDataRows = nan(0,1);
	
    simTimeCol = find(ismember(testVPop.simData.rowInfoNames,'time'));
    interventionIDCol = find(ismember(testVPop.simData.rowInfoNames,'interventionID'));
    elementIDCol = find(ismember(testVPop.simData.rowInfoNames,'elementID'));
    elementTypeCol = find(ismember(testVPop.simData.rowInfoNames,'elementType'));

    % We will also score the VPs to add based on the data
    newVPScores = zeros(1,length(newVPIDs));

    %vpopDataCol1Index = find(ismember(testVPop.expData.Properties.VariableNames,'expVal1'));
    % First column where there is data
    if isa(testVPop,'VPop')
        vpopDataCol1Index = 8;
    else
        vpopDataCol1Index = 12;
    end
    
    addScore = zeros(nTestOutcomes,length(newVPIDs));
    % We'll also calculate a score weighted by biomarker and intervention
    % to avoid biasing towards biomarkers with too many time points
    newVPScoresWBI = zeros(1,length(newVPIDs));
    interventionElement = testVPop.simData.rowInfo(:,[interventionIDCol,elementIDCol,elementTypeCol]);
    [C, ia, interventionElementWeight] = unique(cell2table(interventionElement),'rows');
    % Now create a weight vector
    nInterventionElements = max(interventionElementWeight);
	interventionElementIndices = interventionElementWeight;
    for interventionElementCounter = 1:nInterventionElements
        curRows = find(interventionElementWeight==interventionElementCounter);
        interventionElementWeight(curRows) = 1/length(curRows);
    end
    addScoreWBIWeights = repmat(interventionElementWeight,1,length(newVPIDs));
    
    for rowCounter = 1 : nTestOutcomes
        curSimVals = testVPop.simData.Data(rowCounter,:);
        newVals = curSimVals(newIndices);
        originalValsSort = sort(curSimVals(originalIndices),'ascend');    
        originalValsSort =originalValsSort(find(~isnan(originalValsSort)));
        curSimNonNANIndices = find(~isnan(curSimVals));
        curSimVals = curSimVals(curSimNonNANIndices);
        newNonNANIndices = find(~isnan(newVals));
        newVals = newVals(newNonNANIndices);
        
        simTime = testVPop.simData.rowInfo{rowCounter,simTimeCol};
        interventionID = testVPop.simData.rowInfo{rowCounter,interventionIDCol};
        elementID = testVPop.simData.rowInfo{rowCounter,elementIDCol};
        elementType = testVPop.simData.rowInfo{rowCounter,elementTypeCol};
        % If this is a RECIST VPop, we need to get RECIST-filtered
        % observed experimental data...
        
        expDataRow = find((ismember(testVPop.expData{:,'time'},simTime))&(ismember(testVPop.expData{:,'interventionID'},interventionID))&(ismember(testVPop.expData{:,'elementID'},elementID))&(ismember(testVPop.expData{:,'elementType'},elementType)));
        curExpVals = testVPop.expData{expDataRow,vpopDataCol1Index:end};
        curExpVals = curExpVals(~isnan(curExpVals));

        [~,simValsSortIndices] = sort(curSimVals,'ascend');
        
        originalRange = originalValsSort(length(originalValsSort)) - originalValsSort(1);
        [~, ~, SCall] = alignSamples(sort(curExpVals,'ascend'), sort(curSimVals,'ascend'));

        curBandWidth = (max(SCall) - min(SCall))/10;
        % Calculate VP scores
        % by comparing the cohort PDF to data
        [PDFsim,~] = ksdensity(originalValsSort,SCall,'bandwidth',curBandWidth);%,'Support',[min(SC)-eps max(SC)+eps]);
        [PDFexp,~] = ksdensity(sort(curExpVals,'ascend'),SCall,'bandwidth',curBandWidth);%,'Support',[min(SC)-eps max(SC)+eps]);
        % The PDF integral may not quite be ~1 since the
        % pdf density is not restricted to
        % fall within the observed range
        % We therefore correct the PDFs so we don't bias towards
        % prioritizing the match for biomarkers that have most of their
        % density away from the edges
        PDFsim = PDFsim / trapz(SCall,PDFsim);
        PDFexp = PDFexp / trapz(SCall,PDFexp);
        % We will want to compare the differences 
        pdfDiff = (PDFexp-PDFsim).*(PDFexp>PDFsim);
        
        % These are not needed but included to help debug
        % pdfDiffComb{rowCounter} = pdfDiff;
        % PDFexpComb{rowCounter} = PDFexp;
        % PDFsimComb{rowCounter} = PDFsim;
        % SCallComb{rowCounter}=SCall;
        
        % We want the indices 
        newValInd =nan(1,0);
        scInd =nan(1,0);
        for newCounter = 1 : length(newVals)
            curInd = find(ismember(SCall,newVals(newCounter)));
            if length(curInd) > 0
                scInd = [scInd, curInd(1)];
                newValInd = [newValInd, newCounter];
            end
        end
         
        addScore(rowCounter,newNonNANIndices(newValInd)) = pdfDiff(scInd);
        % We will emphasize VPs that help
        % to address issues with individual
        % variables in the marginal density functions by
        % taking the square
        
    end
    % The sum of the positive difference square terms or just the sum of the
    % positive differences
    % Both have reasonable arguments.
	
    newVPScores = sum(addScore.^2,1);
	% Take the average score for each biomarker/intervention
    newVPScoresWBI = sum(addScore.^2.*addScoreWBIWeights,1);
	% Also take the max score for each biomarker/intervention
	newScoresWBIMax = zeros(nInterventionElements,length(newVPIDs));
	for interventionElementCounter = 1:nInterventionElements
        curRows = find(interventionElementIndices==interventionElementCounter);
		newScoresWBIMax(interventionElementCounter,:) = max(addScore(curRows,:).^2,[],1);
    end 
    newScoresWBIMax = max(newScoresWBIMax,[],1);
    
    if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
        % Reset the RSCORE with PW=1 to get the RSCOREs in the source VPop
        newVPop.pws = (1/length(newVPop.pws)) * ones(1,length(newVPop.pws));
        newVPop = newVPop.addPredTableVals();
        myRTableRECIST = newVPop.rTableRECIST;
        rRowsTarget = newVPop.simData.rRows;
        rRowsSource = find(rRowsTarget>0);         
        
        if ~isempty(rRowsSource)           
            rRowsTarget = rRowsTarget(rRowsSource);      
            [nRRows, nBinCols] = size(myRTableRECIST);          

            curSimValues = nan(nRRows,length(newVPIDs));
            
            myRData = testVPop.simData.rData;
            curSimValues(rRowsTarget,:) = (myRData(rRowsSource, newIndices)); 
            addScore = zeros(nRRows,length(newVPIDs));
            
            for rowCounter = 1 : nRRows
                expPDF = myRTableRECIST{rowCounter,{'expCR','expPR','expSD','expPD'}};
                simPDF = myRTableRECIST{rowCounter,{'predCR','predPR','predSD','predPD'}};
                curWeights = ((expPDF - simPDF) > 0) .* (expPDF - simPDF);
                addScore(rowCounter,:) = (curSimValues(rowCounter,:)==0) .* curWeights(1);
                addScore(rowCounter,:) = addScore(rowCounter,:)+(curSimValues(rowCounter,:)==1) .* curWeights(2);
                addScore(rowCounter,:) = addScore(rowCounter,:)+(curSimValues(rowCounter,:)==2) .* curWeights(3);
                addScore(rowCounter,:) = addScore(rowCounter,:)+(curSimValues(rowCounter,:)==3) .* curWeights(4);
            end
            newVPScoresR = sum(addScore.^2,1);
			%newVPScoresR = max(addScore.^2,[],1);
        end
    end
	
	if flagCheck2D
		[n2DComparisons, ~] = size(newVPop.distTable2D);
		if (n2DComparisons > 0)
			addScore = zeros(n2DComparisons,length(newVPIDs));
            % As a precaution, restart any existing parallel 
            % pools 
%             if ~isempty(gcp('nocreate'))
%                 delete(gcp);
%             end
%             parpool;            
            for rowCounter = 1 : n2DComparisons
                expSample = testVPop.distTable2D{rowCounter,'expSample'}{1};
                simSample = testVPop.distTable2D{rowCounter,'predSample'}{1};
				simPWs = testVPop.distTable2D{rowCounter,'predProbs'}{1};
				% We need to get the predIndices, i.e. the indices that are kept from
				% the original simulation results after applying mechanistic dropouts.
				predIndices = testVPop.distTable2D{rowCounter,'predIndices'}{1};
				[PDFexp, PDFsim, combinedPoints] = align2DPDFs(expSample, simSample, simPWs);
				% We evaluate values from new VPs that haven't dropped out.  i.e.
				% that lie in the intersection of newIndices and predIndices
				newPredIndices = intersect(newIndices,predIndices,'sorted');
                
				% Get the indices for the dropout filtered sample to get the test values
				newPredIndicesFinal = (find(ismember(predIndices,newPredIndices)));
				% simSample is following filtering for dropouts
				testValues = simSample(:,newPredIndicesFinal);
                [~, nTestValues] = size(testValues);
				
				pdfDiff = (PDFexp-PDFsim).*(PDFexp>PDFsim);
				% We want the indices 
				newValInd =nan(1,0);
				scInd =nan(1,0);
                if nTestValues > 0
                    for testCounter = 1 : nTestValues
                        curInd = find(ismember(combinedPoints',(testValues(:,testCounter))','rows'));
                        if length(curInd) > 0
                            % This is the position in combined points
                            scInd = [scInd, curInd(1)];
                            % This is the original index in the new VPs that have not dropped out
                            newValInd = [newValInd, testCounter];
                        end
                    end
                    % Need a map from new VPs that have not dropped out back into new VP position
                    newIndicesNotDropout = find(ismember(newIndices,newPredIndices));

                    % We need to map back to the new VPs and correct
                    % for positions for ones that drop out
                    addScore(rowCounter,newIndicesNotDropout(newValInd)) = pdfDiff(scInd);
                end
            end
%             % Clean up the worker pool
%             delete(gcp)            
			newVPScores2D = sum(addScore.^2,1);
		else
			flagCheck2D = false;
		end
	end
    
	
    if isa(newVPop,'VPopRECIST') || isa(newVPop,'VPopRECISTnoBin')
		if ~flagCheck2D
			newVPScores = [newVPScoresWBI;newVPScoresR;newScoresWBIMax;newVPScores];
		else
			newVPScores = [newVPScoresWBI;newVPScoresR;newVPScores2D;newScoresWBIMax;newVPScores];
		end
    else
        newVPScores = [newVPScoresWBI;newScoresWBIMax;newVPScores];
    end	

    % Select valid VPs from each higher weighted seed VP
    newPassNames = cell(1,0);
    for highCounter = 1 : length(highVPIDs);
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

    % The repeats should already be screened, but just in case
    newPassNames = unique(newPassNames,'first');

    % Finalize the updated cohort
    mergeNames = [originalVPIDs,newPassNames];
    myWorksheet = copyWorksheet(jitteredWorksheet,mergeNames);

    %         % Repopulate results.  Note that this is not fully necessary
    %         % and adds some time to each iteration, but help to avoid issues
    %         % if the original worksheet was carrying old, invalid results.
    %         % This is therefore removed, better to check the worksheet
    %         % VPs before calling expandVPopEffN.
    %         myWorksheet.results = {};
    %         myWorksheet = simulateWorksheet(myWorksheet);
else
    warning(['Unable to proceed in ',mfilename,'.  Returning input worksheet.'])
end
end