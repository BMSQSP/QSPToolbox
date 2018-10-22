function [myWorksheet, newVPop] = expandVPopEffN(myWorksheet,refWorksheet,myExpandVPopEffNOptions,myMapelOptions, hotStartVPop)
% This function takes a worksheet, options object
% mapelOptions object, as well as a reference 
% worksheet for response type bounds,
% tries to fit a VPop. Once successful it adds
% VPs to the worksheet, then tries to fit a VPop with higher effN,
% expanding the effective N in MAPEL in through consecutive iterations.  
% Essentially, it will run MAPEL iteratively, generating VPops until it 
% finds one with an acceptable GOF, generate new VPs based on the most 
% weighted VPs, and select VPs from each most weighted "seed" VP to add.
% Worksheets and VPops are written to file as the algorithm progresses.
% After each iteration, acceptable VPop results and the new
% worksheet are written to file, before the minEffN is increased in MAPEL.
% The algorithm will keep running until it finds a VPop that satisfies the 
% target effN.
% 
% ARGUMENTS:
%  myWorksheet:            A worksheet to expand
%  refWorksheet:           A reference worksheet to use for testBounds.  This should
%                           have the same response types and axes as myWorksheet.
%                           It could be a repeat of myWorksheet if
%                           you just want to maintain current bounds.
%  myExpandVPopEffNOptions An expandVPopEffNOptions object
%  myMapelOptions:         A mapelOptions object
%  hotStartVPop:           A VPop object, will try to hot
%                           start from this solution if provided.
%                           Leave empty or '' if no hot start desired.
%
% RETURNS:
%  myWorksheet
%  newVPop
%
%  NOTES: 
%   The algorithm
%   may encounter problems if an old suffix and
%   wsIterCounter combination are re-used, if existing 
%   VPs in the worksheet have names identical to those 
%   selected by the algorithm.
%

continueFlag = false;
if nargin > 5
    warning(['Too many input arguments to ',mfilename, '. Arguments should be: myWorksheet, refWorksheet, myExpandVPopEffNOptions, and myMapelOptions. Optionally: hotStartVPop.'])
    continueFlag = false;
    hotStartVPop = '';    
elseif nargin > 4
    continueFlag = true;
elseif nargin > 3
    continueFlag = true;
    hotStartVPop = '';
else 
    warning(['Insufficient input arguments to ',mfilename, '. Arguments should be: myWorksheet, refWorksheet, myExpandVPopEffNOptions, and myMapelOptions. Optionally: hotStartVPop.'])
    continueFlag = false;
end

if continueFlag
    if ~strcmp(class(myExpandVPopEffNOptions),'expandVPopEffNOptions')
        warning(['Argument myExpandVPopEffNOptions to ',mfilename, ' should be of class expandVPopEffNOptions.'])
        continueFlag = false;
    end
    if ~ismember(class(myMapelOptions),{'mapelOptions','mapelOptionsRECIST','mapelOptionsRECISTnoBin'})
        warning(['Argument myMapelOptions to ',mfilename, ' should be of class mapelOptions or mapelOptionsRECIST or mapelOptionsRECISTnobin.'])
        continueFlag = false;
    end    
    if ~isequal(myWorksheet.responseTypes,refWorksheet.responseTypes)
        warning(['The response types for myWorksheet and refWorksheet should be the same in ',mfilename, '.'])
        continueFlag = false;
    end         
end 
    



if continueFlag    
    suffix = myExpandVPopEffNOptions.suffix;
    wsIterCounter = myExpandVPopEffNOptions.wsIterCounter;
    targetEffN = myExpandVPopEffNOptions.targetEffN;
    maxNewPerIter = myExpandVPopEffNOptions.maxNewPerIter;
    effNDelta = myExpandVPopEffNOptions.effNDelta;
    minPVal = myExpandVPopEffNOptions.minPVal;
    nTries = myExpandVPopEffNOptions.nTries;
    nRetries = myExpandVPopEffNOptions.nRetries;
    verbose = myExpandVPopEffNOptions.verbose;
    expandCohortSize = myExpandVPopEffNOptions.expandCohortSize;
    restartPVal = myExpandVPopEffNOptions.restartPVal;


    % Starting effN
    curEffN = myMapelOptions.minEffN;

    % Set these for fast optimization
    myMapelOptions.exactFlag = false;
	% Enforce using effN, since
	% we enforce using it for the VPops
    myMapelOptions.useEffN = true;

    % Make sure results are present.  Here, we don't re-simulate VPs if 
    % results are present.
    myWorksheet = simulateWorksheet(myWorksheet);
    refWorksheet = simulateWorksheet(refWorksheet);    
    allResponseTypeIDs = getResponseTypeIDs(myWorksheet);
    nResponseTypes = length(allResponseTypeIDs);

    % Get the responseType evaluation results
    % so we can enforce these as criteria new VPs will need pass.
    % We repopulated these before re-simulating the worksheet.
    myCoeffs = getVPCoeffs(refWorksheet);
    [nAxis, ~] = size(myCoeffs);
    myResponseSummaryTables = cell(1,nResponseTypes);
    testBounds = cell(1,nResponseTypes);
    for responseTypeCounter = 1 : nResponseTypes
        myResponseSummaryTables{responseTypeCounter} = createResponseSummaryTable(refWorksheet, allResponseTypeIDs{responseTypeCounter});
        testBounds{responseTypeCounter} = max(myResponseSummaryTables{responseTypeCounter}.values((nAxis+1):end,:),[],2);
    end

    if (~ismember(class(hotStartVPop),{'VPop','VPopRECIST','VPopRECISTnoBin'}))
        nVPopsFound = 0;
    else
        nVPopsFound = 1;
        oldVPop = hotStartVPop;
        oldVPop.useEffN = true;
        oldVPop.exactFlag = true;
        oldVPop = evaluateGOF(oldVPop);	
    end

    while curEffN <= targetEffN
        bestPVal = 0;
        myMapelOptions.minEffN = curEffN;
        myTestCounter = 0;
        while bestPVal < minPVal
            myTestCounter = myTestCounter + 1;
            myMapelOptions.intSeed=myTestCounter;
            % We will try to restart with initial probabilities, if
            % available, and refresh worksheet data if this is the the first
            % iteration.  We'll also enforce the MAPEL options in
            % case they are tweaked.  Note we don't just directly restart
            % since there may be new simData in the worksheet.
            if (myTestCounter == 1) && (nVPopsFound > 0)
                newVPop = oldVPop;
                newVPop.spreadOut = myMapelOptions.spreadOut;
                newVPop.exactFlag = myMapelOptions.exactFlag;        
                newVPop.useEffN = myMapelOptions.useEffN;    
                newVPop.optimizeTimeLimit = myMapelOptions.optimizeTimeLimit;
                newVPop.optimizeType = myMapelOptions.optimizeType;        
                newVPop.optimizePopSize = myMapelOptions.optimizePopSize;
                newVPop.intSeed = myMapelOptions.intSeed;
                newVPop.nIters = myMapelOptions.nIters;
                newVPop.tol = myMapelOptions.tol;
                newVPop.expData = myMapelOptions.expData;
                newVPop.mnSDTable = myMapelOptions.mnSDTable;
                newVPop.binTable = myMapelOptions.binTable;
                newVPop.distTable = myMapelOptions.distTable; 
                if isa(myMapelOptions,'mapelOptionsRECIST') || isa(myMapelOptions,'mapelOptionsRECISTnoBin')
                    newVPop.distTable2D = myMapelOptions.distTable2D;
                    newVPop.relSLDvar = myMapelOptions.relSLDvar;
                    newVPop.absALDVar = myMapelOptions.absALDVar;
                    newVPop.crCutoff = myMapelOptions.crCutoff;                      
                    newVPop.brTableRECIST = myMapelOptions.brTableRECIST;
                    newVPop.rTableRECIST = myMapelOptions.rTableRECIST;                    
                    newVPop.recistSimFilter = createRECISTSimFilter(myWorksheet, newVPop);                    
                end		
				newVPop.minEffN = myMapelOptions.minEffN;
				newVPop = newVPop.getSimData(myWorksheet);
				if ~isa(newVPop,'VPopRECISTnoBin')
					newVPop = newVPop.assignIndices(myWorksheet, myMapelOptions);
					% We assume the old bin probs are compatible
					% with the new indices assignments.  That is,
					% there are an equal number of bins and axes.
					newVPop.binProbs = oldVPop.binProbs;
				else
					% Otherwise we are weighting on VPs directly.
					% We need to assign pws for added vps, also
					% the CoeffsTable
					newVPop=newVPop.assignCoeffs(myWorksheet);
					[~,nVP_nobin]=size(newVPop.coeffsTable);
					nVP_diff=nVP_nobin-length(oldVPop.pws);
					if nVP_diff > 0
						newIndices = [length(oldVPop.pws)+1 : nVP_nobin];
						allVPIDs = getVPIDs(myWorksheet);
						parentIndices = nan(1,nVP_diff);
						for index = 1 : nVP_diff
							curVPID = allVPIDs{newIndices(index)};
							parentID = strsplit(curVPID,['_']);
							nIndices = length(parentID);
							parentID = parentID(1:nIndices-2);
							parentID = strjoin(parentID,'_');
							parentIndices(index) = find(ismember(allVPIDs,parentID));
							% This will error out if the parent index is > 1.
							% That should not happen!
						end
						uniqueIndices = unique(parentIndices);
						newPWs = nan(1,nVP_diff);
                        oldPWs = oldVPop.pws;
						for parentCounter = 1 : length(uniqueIndices)
							parentIndex = uniqueIndices(parentCounter);
							totalWeight = oldVPop.pws(parentIndex);
							childrenIndices = find(ismember(parentIndices,parentIndex));
							nChildren = length(childrenIndices);
							oldPWs(parentIndex) = totalWeight/(nChildren+1);
							newPWs(childrenIndices) = totalWeight/(nChildren+1);
						end
						newVPop.pws=[oldPWs, newPWs];
					else
						pws_addition=sum(ones(1,nVP_diff)*1/nVP_nobin);
						newVPop.pws=[oldVPop.pws*(1-pws_addition), ones(1,nVP_diff)*1/nVP_nobin];  	
					end
					% **method 1: set addtional vps has 0 pws
	%                 newVPop.pws=[oldVPop.pws,zeros(1,nVP_diff)]; % maybe not so good?
					% ** end of method 1
					% ++method 2: divide the highest pws and distribute
	%                 [pws_max_temp, pws_max_temp_ind]=max(oldVPop.pws);
	%                 pws_max_avg=pws_max_temp/(nVP_diff+1);
	%                 oldVPop.pws(pws_max_temp_ind)=pws_max_avg;
	%                 newVPop.pws=[oldVPop.pws,ones(1,nVP_diff)*pws_max_avg];
					% ++ end of method 2
					% method 3: assign 1/nVP_nobin as the weights of new VPs
					% and rescale all others
					% method 4: distribute the parent weight to the children
                end
                % We randomize here with magnitude myExpandVPopEffNOptions.expandRandomStart
				% If the magnitude is set to infinity, we just fully randomize before "restarting"
				if isinf(myExpandVPopEffNOptions.expandRandomStart)
					if isa(newVPop, 'VPopRECISTnoBin')
						newVPop = newVPop.startPWs(myWorksheet,true);
					else
						newVPop = newVPop.startProbs(true)
					end	
					newVPop = restartMapel(newVPop, 0);
				else
					newVPop = restartMapel(newVPop, myExpandVPopEffNOptions.expandRandomStart);
				end
            else
				% For consistency, on the first iteration
                % if a VPop isn't provided we will use the
				% myExpandVPopEffNOptions.expandRandomStart
				% value to decide randomization status.
				curMapelOptions = myMapelOptions;
				if (myTestCounter == 1) 
					curMapelOptions.randomStart = myExpandVPopEffNOptions.expandRandomStart;
				end
                newVPop = mapel(myWorksheet, curMapelOptions);
            end
            newVPop.useEffN = true;
            newVPop.exactFlag = true;
            newVPop = evaluateGOF(newVPop);
            curVPopEffN = 1/sum(newVPop.pws.^2);
            if verbose
                disp(['Previous best pvalue ',num2str(bestPVal),' with effN ',num2str(curEffN),'.'])
                disp(['Current solution effN ', num2str(curVPopEffN), ' with pvalue ',num2str(newVPop.gof),' in ',mfilename,'.'])
                disp(['Current iteration ', num2str(myTestCounter),'.'])
            end
            if ((newVPop.gof < minPVal) && (newVPop.gof >= restartPVal) && (curVPopEffN >= curEffN))
                disp(['This pvalue not quite at acceptance, restarting MAPEL with this as initial guess.  Using myMapelOptions.randomStart value to set the magnitude of randomization.'])
                myRandomStart = myMapelOptions.randomStart;
                restartCounter = 1;
                if (newVPop.gof > bestPVal)
					% Save here in case this is the best
                    mySaveName = ['vpop_',suffix,'_iter',num2str(wsIterCounter),'_effN',num2str(curEffN),'_curBest'];
                    saveVPop(newVPop, mySaveName);
                    bestPVal = newVPop.gof;   
                end								
                while (((restartCounter <= nRetries) && ~(newVPop.gof >= minPVal)) || ((newVPop.gof > (0.95 * minPVal)) && ~(newVPop.gof >= minPVal)))
                    restartVPop = newVPop;
                    restartVPop.useEffN = myMapelOptions.useEffN;
                    restartVPop.exactFlag = myMapelOptions.exactFlag;
                    restartVPop.intSeed = myTestCounter + restartCounter;
                    % Increase the tolerance stringency for restarting
                    restartVPop.tol = newVPop.tol*0.1;
                    restartVPop = restartMapel(restartVPop, myRandomStart);
                    restartVPop.tol = newVPop.tol;
                    restartVPop.useEffN = true;
                    restartVPop.exactFlag = true;
                    restartVPop = evaluateGOF(restartVPop);				
                    restartVPopEffN = 1/sum(restartVPop.pws.^2);
                    if verbose
                        disp(['Previous best pvalue ',num2str(bestPVal),' with effN ',num2str(curEffN),'.'])
                        disp(['Current restart solution effN ', num2str(restartVPopEffN), ' with pvalue ',num2str(restartVPop.gof),' in ',mfilename,'.'])
                        disp(['Current restart iteration ', num2str(restartCounter),'.'])
                    end
                    restartCounter = restartCounter + 1;
                    if ((restartVPop.gof >= newVPop.gof) && (restartVPopEffN >= curEffN))
                        newVPop = restartVPop;
                        curVPopEffN = restartVPopEffN;
                        if (newVPop.gof > bestPVal)
                            % We also do this during the restart loop so we can
                            % track progress
                            mySaveName = ['vpop_',suffix,'_iter',num2str(wsIterCounter),'_effN',num2str(curEffN),'_curBest'];
                            saveVPop(newVPop, mySaveName);
                            bestPVal = newVPop.gof;   
                        end
                    end
                end
            end
            if ((newVPop.gof > bestPVal) && (curVPopEffN >= curEffN))
                mySaveName = ['vpop_',suffix,'_iter',num2str(wsIterCounter),'_effN',num2str(curEffN),'_curBest'];
                saveVPop(newVPop, mySaveName);
                bestPVal = newVPop.gof;
            end
            % If we have surpassed nTries VPop iterations without success 
            % we will try adding VPs based on the last, best VPop.
            if (mod(myTestCounter,nTries) == 0) && (nVPopsFound > 0) && ~((newVPop.gof > minPVal) && (1/sum(newVPop.pws.^2) >= curEffN))
                if expandCohortSize > 0
                    wsIterCounter = wsIterCounter + 1;
                    [myWorksheet, newPassNames] = expandWorksheetVPsFromVPop(myWorksheet,oldVPop, myMapelOptions,suffix,wsIterCounter, maxNewPerIter, testBounds, expandCohortSize, myExpandVPopEffNOptions.gaussianStd, myExpandVPopEffNOptions.maxNewPerOld);
                    saveWorksheet(myWorksheet,['myWorksheet_',suffix,'_iter',num2str(wsIterCounter)]); 
                    if verbose
                        disp(['Unable to find an acceptable VPop with initial worksheet in ',num2str(nTries),' VPop fit restarts, added ', num2str(length(newPassNames)), ' VPs to the worksheet in ',mfilename,' to start worksheet iteration ',num2str(wsIterCounter),'.'])
                    end       
                end                
            end
        end
        nVPopsFound = nVPopsFound+1;
        oldVPop = newVPop;
        % Increment target
        curEffN = curEffN + effNDelta;
        myMapelOptions.minEffN = curEffN;

        if curEffN <= targetEffN
            if expandCohortSize > 0
                % Get ready to expand VPs
                wsIterCounter = wsIterCounter+1;
                [myWorksheet, newPassNames] = expandWorksheetVPsFromVPop(myWorksheet,oldVPop, myMapelOptions,suffix,wsIterCounter, maxNewPerIter, testBounds, expandCohortSize, myExpandVPopEffNOptions.gaussianStd, myExpandVPopEffNOptions.maxNewPerOld);
                % Save the new worksheet that will be used.
                saveWorksheet(myWorksheet,['myWorksheet_',suffix,'_iter',num2str(wsIterCounter)]);   
                if verbose
                    disp(['Added ', num2str(length(newPassNames)), ' VPs to the worksheet in ',mfilename,' to start worksheet iteration ',num2str(wsIterCounter),'.'])
                end
            end
        end
    end
else
    warning(['Exiting ',mfilename, '. Returning input worksheet and hot start VPop.'])    
    newVPop = hotStartVPop;
end
end