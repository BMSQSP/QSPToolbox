function myVPop = mapel(myWorksheet, myMapelOptions)
% This is the main function to call MAPEL, which will develop a virtual
% population based on the simulation results in the worksheet and the
% settings in the mapelOptions.
%
% ARGUMENTS
%  myWorksheet:            (Required) a worksheet instance.  We allow a 
%                          special pass-through where a VPop object can also
%                          be given here to get a couple critical
%                          properties (simData, indexTable, binEdges
%                          binMidPoints) mapel to be called from
%                          a restart, but in this case all other
%                          properties are taken from the options.
%  myMapelOptions:         (Required) an instance of a mapelOptions object.
%
% % Returns
%  VPop:                    an instance of a VPop
%

continueFlag = false;
if nargin > 2
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: myWorksheet and a mapelOptions (or mapelOptionsRECIST or mapelOptionsRECISTnoBin) object'])
    continueFlag = false;
elseif nargin > 1
    continueFlag = true;
else
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: myWorksheet and a mapelOptions (or mapelOptionsRECIST or mapelOptionsRECISTnoBin) object.'])
end

if continueFlag
    if ~ismember(class(myMapelOptions),{'mapelOptions','mapelOptionsRECIST','mapelOptionsRECISTnoBin'})
        continueFlag = false;
        warning(['Input mapelOptions not recognized in call to ',mfilename,'.  Requires: myWorksheet and a mapelOptions or mapelOptionsRECIST or mapelOptionsRECISTnoBin object.'])
    end       
end

if continueFlag
    % First copy all of the needed properties over
    myVPop = initializeOptionPropertiesToVPop(myMapelOptions);
    
    % Now add additional properties not in the mapelOptions
    if isa(myVPop,'VPopRECIST')      
        if ~isa(myWorksheet,'VPopRECIST')
            myVPop.recistSimFilter = createRECISTSimFilter(myWorksheet, myVPop);
        else   
            myVPop.recistSimFilter = myWorksheet.recistSimFilter;
        end
    end	    
	if isa(myVPop,'VPopRECISTnoBin')      
        if ~isa(myWorksheet,'VPopRECISTnoBin')
            myVPop.recistSimFilter = createRECISTSimFilter(myWorksheet, myVPop);
        else   
            myVPop.recistSimFilter = myWorksheet.recistSimFilter;
        end
    end		 

    % Next, we create a table of nVPs x nBins to index each VP axis into a
    % bin for subsequent calculations.
    % We go to the simulation results in the worksheet and create
    % a table of the needed results to compare to the experimental data
    if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPop')
        if ~isa(myWorksheet,'VPop') && ~isa(myWorksheet,'VPopRECIST')
            myVPop = myVPop.assignIndices(myWorksheet, myMapelOptions);
            myVPop = myVPop.getSimData(myWorksheet);
			myVPop = myVPop.assignCoeffs(myWorksheet);
            if myVPop.spreadOut > 0
                myVPop.coeffsDist = single(pdist2(myVPop.coeffsTable',myVPop.coeffsTable'));
            end
        else
            myVPop.simData = myWorksheet.simData;
            myVPop.indexTable = myWorksheet.indexTable;
            myVPop.binEdges = myWorksheet.binEdges;
            myVPop.binMidPoints = myWorksheet.binMidPoints;
			myVPop.coeffsTable = myWorksheet.coeffsTable;
            if myVPop.spreadOut > 0
                myVPop.coeffsDist = single(pdist2(myVPop.coeffsTable',myVPop.coeffsTable'));
            end
        end
    else %'VPopRECISTnoBin'
        if ~isa(myWorksheet,'VPopRECISTnoBin')
            myVPop = myVPop.getSimData(myWorksheet);
            myVPop = myVPop.assignCoeffs(myWorksheet);
            if myVPop.spreadOut > 0
                myVPop.coeffsDist = single(pdist2(myVPop.coeffsTable',myVPop.coeffsTable'));
            end
        else
            myVPop.simData = myWorksheet.simData;
            % PWs will still be taken from the mapelOptionsNoBins
            myVPop.coeffsTable = myWorksheet.coeffsTable;
            if myVPop.spreadOut > 0
                myVPop.coeffsDist = single(pdist2(myVPop.coeffsTable',myVPop.coeffsTable'));
            end
        end
    end
    												
    if myVPop.intSeed > -1
        rng(myVPop.intSeed, 'twister');
    end        
	
    % Also assign the now fixed sim data into the bin tables initially
    % rather than during execution to reduce table assignments and
    % recalculating the 1D mesh.
    myVPop = myVPop.addDistTableSimVals();
    
    myVPop = myVPop.addDistTable2DSimVals();
	myVPop = myVPop.addCorTableSimVals();
    	
    
    if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPop')
        % We adopt the index table convention from the original MAPEL paper
        myIndexTable = myVPop.indexTable;
        [myNAxis, nVP] = size(myIndexTable);
        myBinEdges = myVPop.binEdges;
        myBinMidPoints = myVPop.binMidPoints;
        myNBins = myMapelOptions.nBins;

        % We can continue a run from a previous run if we have a valid
        % initial probability table.
        myInitialProbs = myMapelOptions.initialProbs;
        myRandomStart = myMapelOptions.randomStart;

        if isequal([myNAxis,myNBins], size(myInitialProbs))
            if myRandomStart > 0
                for axisCounter = 1 : myNAxis
                    myInitialProbsTrans = hyperTransform(myInitialProbs(axisCounter,:));
                    mySDs = myRandomStart*abs(myInitialProbsTrans);
                    myInitialProbsTrans = myInitialProbsTrans + randn(1,myNBins-1).*mySDs;
                    myInitialProbs(axisCounter,:) = invHyperTransform(myInitialProbsTrans); 
                end
            end
            myVPop.binProbs = myInitialProbs;
            myVPop = myVPop.assignPWs();
        else
            myVPop = myVPop.startProbs(myRandomStart>0);  
            myVPop = myVPop.assignPWs();
        end
    else 
		[nAxis, nVP] = size(myVPop.coeffsTable);
        
		myInitialPWs = myMapelOptions.initialPWs;
		myRandomStart = myMapelOptions.randomStart;
		
		if isequal(1, length(myInitialPWs))
			if myInitialPWs < 0
				% This implies we want to try starting from near an
				% optimal solution to the linearized problem
				myOptimOptions = LinearCalibrationOptions();
				myOptimOptions.cdfProbsToFit = 0.05:0.05:0.95;
				myOptimOptions.optimizationAlgorithm = "nnls";
				myOptimOptions.optimizationAlgorithmOptions.Accy = 0;
				myOptimOptions.priorPrevalenceWeightAssumption = 'specified';
                % Need to add predTableVals to run.
                myVPop = myVPop.startPWs(myWorksheet,0);
                myVPop = myVPop.addPredTableVals();
				
				linearCalibrationObject = LinearCalibration(myVPop,'optimOptions',myOptimOptions);
				linearCalibrationObject = linearCalibrationObject.run();
				
				if linearCalibrationObject.OptimizationResults.exitFlag == 1
                    % Run 2x as results are sensitive to prior PW
                    % assumption
					linearCalibrationObject = LinearCalibration(linearCalibrationObject.OptimizedVPop,'optimOptions',myOptimOptions);
					linearCalibrationObject = linearCalibrationObject.run();					
					myInitialPWs = linearCalibrationObject.OptimizedVPop.pws;
					if linearCalibrationObject.OptimizationResults.exitFlag == 1
						myInitialPWs = linearCalibrationObject.OptimizedVPop.pws;
					end
				else
					warning(['Unable to find optimal solution to linear problem in ',mfilename,'.  Proceeding with default options.'])
					myInitialPWs = [];
				end
			end
		end
        
		if isequal(nVP, length(myInitialPWs))
			if myRandomStart > 0
				myInitialProbsTrans = hyperTransform(myInitialPWs);
				mySDs = myRandomStart*abs(myInitialProbsTrans);
				myInitialProbsTrans = myInitialProbsTrans + randn(1,nVP-1).*mySDs;
				myInitialPWs = invHyperTransform(myInitialProbsTrans); 
			end
			myVPop.pws = myInitialPWs;
		else
			myVPop = myVPop.startPWs(myWorksheet,myRandomStart>0);   
		end
		
    end    

    % Update table values.  Not strictly necessary but a nice
    % step to diagnose and does have much computational cost            
    myVPop = myVPop.addPredTableVals();    
    myVPop = findFit(myVPop);
else
    myVPop = VPop;
    warning(['Unable to run ',mfilename,'.  Returning an empty VPop.'])
end
end