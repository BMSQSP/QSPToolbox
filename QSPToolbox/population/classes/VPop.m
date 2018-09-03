classdef VPop
% This is a virtual population class.  This contains the results if a call
% to MAPEL.  There should be too much here to that needs to
% be changed directly, with the biggest possible exception being adjusting 
% useEffN before re-evaluating the GOF.  However, there is a substantial
% amound of customization possible by resuming fits with an existing
% VPop so you may find reasons to modify the properties.
%
% PROPERTIES:
%  indexTable:   (Don't manipulate) A nAxis X nVP table that indicates in which 
%                axis bin each VP falls into.  This is usually populated at 
%                the beginning of the MAPEL algorithm by the assignIndices
%                method.
%  binEdges:     (Don't manipulate) Edges for the axis bins. These are also
%                populated by the assignIndices method.
%  binMidPoints: (Don't manipulate) Mid points for the axis bins, populated
%                by the assignIndices method.
%  binProbs:     (Don't manipulate) Probabilities for each axis bin.  These are directly
%                varied in the optimization, and the results are used
%                to calculate the individual prevalence weights.
%  pws:          (Don't manipulate) A 1xnVP vector of prevalence weights,
%                which are calculated based on the assignPWs.
%  expData:      (Don't manipulate) A table of experimental data used to
%                guide optimization.  It is usually taken from the
%                mapelOptions.
%  simData:      (Don't manipulate) A table of simulation results that 
%                pairs to the experimental data.  Usually assigned in
%                MAPEL.
%  mnSDTable:    (Don't manipulate) Summary table of experimental and 
%                simulated (weighted) mean and standard deviation.  Usually
%                copied from the mapelOptions and updated with
%                the weighted simulation results from the VPop fitting.
%                If you want to adjust a fit and restart, you may want to 
%                adjust some of the weights here.
%  binTable:     (Don't manipulate) Summary table of experimental and  
%                simulated binned (weighted) data for a virtual population
%                Usually copied from the mapelOptions and updated with
%                the weighted simulation results from the VPop fitting.
%                If you want to adjust a fit and restart, you may want to 
%                adjust some of the weights here.
%  distTable:    (Don't manipulate) Summary table of experimental and  
%                simulated distribution data for a virtual population
%                Usually copied from the mapelOptions and updated with
%                the weighted simulation results from the VPop fitting.
%                If you want to adjust a fit and restart, you may want to 
%                adjust some of the weights here. 
%  gofMn:        (Don't manipulate) Goodness of fit statistics for 
%                individual endpoint/time mean comparisons between data
%                and virtual population.  Usually calcualated by the
%                evaluateGOF function.
%  gofSD:        (Don't manipulate) Goodness of fit statistics for 
%                individual endpoint/time standard deviation comparisons  
%                between data and virtual population.  Usually calculated 
%                by the evaluateGOF function.
%  gofBin:       (Don't manipulate) Goodness of fit statistics for 
%                individual endpoint/time bin distribution comparisons 
%                between data and virtual population.  Usually calculated 
%                by the evaluateGOF function.
%  gofDist:      (Don't manipulate) Goodness of fit statistics for 
%                empirical distribution comparisons 
%                between data and virtual population.  Usually calculated 
%                by the evaluateGOF function.
%  gof:          (Don't manipulate) Composite goodness-of-fit result.
%                Usually updated by a call to the evaluateGOF.
%  useEffN:      Whether or not to use the effN for optimization and
%                evaluation of the GOF.  It is generally recommended to
%                useEffN for evaluation of the final fit but not
%                during optimization.
%  exactFlag:    Whether to check to apply Fisher's exact test instead
%                of the chi-square approximation for binned distribution
%                test.  This may be strictly correct but can slow
%                calculations, often with small impacts on GOF
%                calculations.
%  spreadOut:    Used for optimization, if nonzero, will  
%                apply a penalty proportional to spreadOut to increase how
%                the prevalence weight is distributed.
%  optimizeTimeLimit:   Time limit for optimizing the VPop in s.
%  optimizeType:        Type of optimization algorithm to employ: "pso,"
%                       "ga," or "simplex".  Default is
%                       "pso".
%  optimizePopSize:     Number of solutions to try in each optimization
%                       generation.  Only applies to "pso" and "ga".  This 
%                       is the population size to use in the optimization 
%                       runs. Default is 1000.
%  intSeed:             A non-negative integer seed to initialize the 
%                       random number generator.  Set to -1 to avoid
%                       changing the state of the random number generator.
%  tol:                 Numerical tolerance for optimization.
%  nIters:              Maximum number of iterations for fminsearch.
%  minEffN:             Minimum effective N.  A large penalty is applied
%                       if the effN drops below this N during optimization,
%                       better ensuring solutions that weight multiple VPs.
%                       The minEffN setting operates independently of 
%                       useEffN, and often works better than modifying 
%                       spreadOut.  The default is 0.
%

properties
    indexTable
    binEdges
    binMidPoints
    binProbs
    pws
    expData
    simData      
    % We combine simulation and experimental results
    % into the same table.
    mnSDTable
    binTable
    distTable
    gofMn
    gofSD
    gofBin
    gofDist
    gof
    useEffN  
    exactFlag
    spreadOut     
    optimizeTimeLimit
    optimizeType
    optimizePopSize
    intSeed  
    tol
    nIters
    minEffN
end

methods
      function obj = set.indexTable(obj,myIndexTable)
          obj.indexTable = myIndexTable;
      end

      function obj = set.binEdges(obj,myBins)
          obj.binEdges = myBins;
      end

      function obj = set.binMidPoints(obj,myBinMidPoints)
          obj.binMidPoints = myBinMidPoints;
      end

      function obj = set.binProbs(obj,myBinProbs)
          % We enforce constraints for binProbs:
          % a minimum value is 1E-14 rows must sum to 1
          % This helps avoid issues with the solver as well
          [nrows,ncols] = size(myBinProbs);
          if ((nrows > 0) && (ncols > 0))
            myBinProbs = myBinProbs ./ (repmat(sum(myBinProbs,2),1,ncols));
            adjustFlag = myBinProbs < 1E-14;
            myNadjust = sum(adjustFlag,2);
            myValDecrease = 1E-14 ./ myNadjust; 
            myValDecrease(isinf(myValDecrease)) = 0 ;
            myBinProbs = (myBinProbs + adjustFlag * 1E-14) - ((1- adjustFlag) .* (repmat(myValDecrease,1,ncols)));
            myBinProbs = myBinProbs ./ (repmat(sum(myBinProbs,2),1,ncols));
          end
          obj.binProbs = myBinProbs;
      end

      function obj = set.pws(obj,myPWs)
          obj.pws = myPWs;
      end

      function obj = set.mnSDTable(obj,myMnSDTable)
          obj.mnSDTable = myMnSDTable;
      end

      function obj = set.binTable(obj,myBinTable)
          obj.binTable = myBinTable;
      end

      function obj = set.distTable(obj,myDistTable)
          obj.distTable = myDistTable;
      end      
      
      function obj = set.expData(obj,myExpData)
          obj.expData = myExpData;
      end

      function obj = set.simData(obj,mySimData)
          obj.simData = mySimData;
      end

      function obj = set.gofMn(obj,myGOF)
          obj.gofMn = myGOF;
      end

      function obj = set.gofSD(obj,myGOF)
          obj.gofSD = myGOF;
      end
      
      function obj = set.gofBin(obj,myGOF)
          obj.gofBin = myGOF;
      end     
      
      function obj = set.gofDist(obj,myGOF)
          obj.gofDist = myGOF;
      end          
      
      function obj = set.gof(obj,myGOF)
          obj.gof = myGOF;
      end          
      
      function obj = set.useEffN(obj,myUseEffN)
          if islogical(myUseEffN)
               obj.useEffN = myUseEffN;
          else
               error(['Property useEffN in ',milename,' must be logical.'])
          end
      end  
      
      function obj = set.exactFlag(obj,myFlag)
          obj.exactFlag = myFlag;
      end       
      
      function obj = set.spreadOut(obj,mySpreadOut)
          if (mySpreadOut >= 0)
               obj.spreadOut = mySpreadOut;
          else
               error(['Property spreadOut in ',milename,' must be >= 0.'])
          end
      end    
      
      function obj = set.optimizeTimeLimit(obj,myOptimizeTimeLimit)
          obj.optimizeTimeLimit = myOptimizeTimeLimit;
      end  
      
      function obj = set.optimizeType(obj,myOptimizeType)
          if sum(ismember({'pso','ga','simplex'},lower(myOptimizeType))) == 1
              obj.optimizeType = lower(myOptimizeType);
          else
              error(['Property optimizeType in ',milename,' must be "ga," "pso," or "simplex"'])
          end
      end          
      
      function obj = set.optimizePopSize(obj,myPopulationSize)
          if ((isnumeric(myPopulationSize) == true) && (myPopulationSize >= 0))
              obj.optimizePopSize = myPopulationSize;
          else
            error('Invalid optimizePopSize value')
          end
      end             
      
      function obj = set.intSeed(obj,myIntSeed)
          if ((mod(myIntSeed,1) == 0) && (myIntSeed>=-1))
              obj.intSeed = (myIntSeed);
          else
              error(['Invalid intSeed specified for ',mfilename,', a non-negative integer should be specified, or -1 to ignore.'])
          end
      end      
      
      function obj = set.tol(obj,myTol)
          if (myTol > 0)
               obj.tol = myTol;
          else
               error(['Property tol in ',milename,' must be > 0.'])
          end
      end       
      
      function obj = set.nIters(obj,myNIters)
          if (myNIters > 0) && (mod(myNIters,1) == 0)
               obj.nIters = myNIters;
          else
               error(['Property nIters in ',milename,' must be a positive integer.'])
          end
      end         
      
      function obj = set.minEffN(obj,myMinEffN)
          if ((isnumeric(myMinEffN) == true) && (myMinEffN >= 0))
              obj.minEffN = myMinEffN;
          else
            error('Invalid minEffN value')
          end
      end       
      
      function value = get(obj,propName)
          switch propName
              case 'indexTable'
                  value = obj.indexTable;
              case 'binEdges'
                  value = obj.binEdges;
              case 'binMidPoints'
                  value = obj.binMidPoints;
              case 'binProbs'
                  value = obj.binProbs;   
              case 'pws'
                  value = obj.pws; 
              case 'mnSDTable'
                  value = obj.mnSDTable; 
              case 'binTable'
                  value = obj.binTable; 
              case 'distTable'
                  value = obj.distTable;                   
              case 'expData'
                  value = obj.expData;                   
              case 'simData'
                  value = obj.simData;            
              case 'gofMn'
                  value = obj.gofMn;      
              case 'gofSD'
                  value = obj.gofSD;   
              case 'gofBin'
                  value = obj.gofBin;   
              case 'gofDist'
                  value = obj.gofDist;                     
              case 'gof'
                  value = obj.gof;   
              case 'spreadOut'
                  value = obj.spreadOut;
              case 'useEffN'
                  value = obj.useEffN;
              case 'exactFlag'
                  value = obj.exactFlag;                   
              case 'optimizeTimeLimit'
                  value = obj.optimizeTimeLimit;
              case 'optimizeType'
                  value = obj.optimizeType;  
              case 'optimizePopSize'
                  value = obj.optimizePopSize;
              case 'intSeed'
                  value = obj.intSeed;             
              case 'tol'
                  value = obj.tol;
              case 'nIters'
                  value = obj.nIters;       
              case 'minEffN'
                  value = obj.minEffN;                    
              otherwise
                  error(['Error: ',propName ,' is not a valid ',mfilename,' property.'])
          end
      end

      function obj = assignIndices(obj, myWorksheet, myMapelOptions)
          % This method takes a worksheet and mapelOptions structure and
          % creates a bin table and other poperties needed for the virtual
          % population and MAPEL algorithm.
          %
          % ARGUMENT
          %  (self)
          %  myWorksheet:    a worksheet structure, with VPs, axes, and
          %                  coefficients
          %  myMapelOptions: an MAPEL options object instance, with
          %                  properties
          %                   nBins
          %                   equalBinBreaks
          %
          % RETURNS
          %  (self): properties of the virtual population are updated:
          %           indexTable
          %           binEdges
          %           binMidPoints
          %  
          
          myVPCoeffs = getVPCoeffs(myWorksheet);
          [nAxis, nVP] = size(myVPCoeffs);
          myIndexTable = ones(nAxis, nVP);
          nBins = myMapelOptions.nBins;
          myBinEdges = nan(nAxis, (nBins+1));
          myBinMidPoints = nan(nAxis, nBins);
          equalBinBreaks = myMapelOptions.equalBinBreaks;
          for axisCounter = 1 : nAxis
              % We currently only implement continuous variables from
              % the paper/R MAPEL   
              if equalBinBreaks
                   % Rather than scaling based on the sampled min and
                   % max as in the R MAPEL algorithm, we adjust according
                   % to the allowed,
                   % and this range is [0, 1] by the axis definition.
                   myBinEdges(axisCounter,:) =  (0:1/nBins:1);
              else
                   myPercentiles =  (0:1/nBins:1)';
                   myBinEdges(axisCounter,:) = (quantile(myVPCoeffs(axisCounter,:),myPercentiles));
              end
              % We'll use the cdf convention that FX(x) = P(X <= x), so 
              % check if a value is <= the bin upper cutoff in order to 
              % assign it.
              % Note that effectively we want to ignore the first bin of 0
              % as a cutoff and lump this in with the first bin.
              for binCounter = 1 : nBins
                   myBinMidPoints(axisCounter,binCounter) = (myBinEdges(axisCounter,binCounter) + myBinEdges(axisCounter,binCounter+1))/2;
                   % First bins are assigned an index of 0 by this method.
                   myIndexTable(axisCounter, :) = myIndexTable(axisCounter, :) + (myVPCoeffs(axisCounter, :) > myBinEdges(axisCounter,binCounter+1));
              end     
          end
          obj.indexTable = myIndexTable;
          obj.binEdges = myBinEdges;
          obj.binMidPoints = myBinMidPoints;
      end

      function obj = startProbs(obj, myRandomStart)
          % This method initializes bin probabilities, and is used
          % prior to the optization.
          %
          % ARGUMENTS
          %  (self)
          %  myRandomStart:  A boolean variable (true/false), if true
          %                  the bin probabilities will be set
          %                  randomly and if false they will be uniform
           if nargin < 1
                myRandomStart = false;
           end
           myBinMidPoints = obj.binMidPoints;
           [nAxes, nBins] = size(myBinMidPoints);
           if ~myRandomStart
                myUniformStartProbs = ones(nAxes, nBins) ./ nBins;
           else
                myUniformStartProbs = rand([nAxes, nBins]);
                for axisCounter = 1 : nAxes
                    myUniformStartProbs(axisCounter,:) = myUniformStartProbs(axisCounter,:)/sum(myUniformStartProbs(axisCounter,:));
                end
           end
           obj.binProbs = myUniformStartProbs;
      end

      function obj = assignPWs(obj)
          % Update prevalence weight assignment for the virtual population.
          % Note that existing predicted population results will be set
          % to NaN in:
          %  mnSDTable
          %  binTable
          %  distTable
          %
          % And also the following properties will be reset:
          %  gofMn
          %  gofSD
          %  gofBin
          %  gofDist
          %  gof
          %
          % ARGUMENTS:
          %  (self): No additional arguments, but the VPop object must have
          %          previously assigned the following properties:
          %           indexTable
          %           binProbs
          %
          % RETURNS:
          %  (self): Returns the VPop with the updated PW properties.
          %
           myIndexTable = obj.indexTable;
           [myNAxis, myNVP] = size(myIndexTable);
           myPWs = zeros(1,myNVP);
           myBinProbs = obj.binProbs;
           %pwCutoff = obj.pwCutoff;
           for axisCounter = 1 : myNAxis
               % Note this will result in -Inf if any bins have
               % prob of zero.  MATLAB's exponential at the end corrects
               % this, but there is a danger of getting back all
               % -inf if all of the bin probabilities are small.
                myPWs = myPWs + log(myBinProbs(axisCounter,myIndexTable(axisCounter,:)));
           end
           % Carried over from paper MAPEL: force
           % biggest value to ~1 to avoid underflow
           myPWs = myPWs - max(myPWs);
           myPWs = exp(myPWs); %/ sum(exp(myPWs));
           %myPWs = myPWs.*(myPWs >= pwCutoff);
           myPWs = myPWs / sum(myPWs);
           obj.pws = myPWs;
           
           % If we re-calculate PWs, eliminate any
           % previous PW-dependent properties - i.e. those
           % that start with "pred"
           myMnSDTable = obj.mnSDTable;
           [nRows, nCols] = size(myMnSDTable);
           if nRows > 0
                myString = 'pred';
                varNames = myMnSDTable.Properties.VariableNames;
                myPos = find(strncmpi(myString,varNames,length(myString)));
                nanifyVars = varNames(myPos);
                for varCounter = 1 : length(nanifyVars)
                    myMnSDTable.(nanifyVars{varCounter}) = nan(nRows,1);
                end
           end
           obj.mnSDTable = myMnSDTable;

           myBinTable = obj.binTable;
           [nRows, nCols] = size(myBinTable);
           if nRows > 0
                myString = 'pred';
                varNames = myBinTable.Properties.VariableNames;
                myPos = find(strncmpi(myString,varNames,length(myString)));
                nanifyVars = varNames(myPos);
                for varCounter = 1 : length(nanifyVars)
                    myBinTable.(nanifyVars{varCounter}) = nan(nRows,1);
                end
           end
           obj.binTable = myBinTable;
           
           myDistTable = obj.distTable;
           [nRows, nCols] = size(myDistTable);
           if nRows > 0
                myNanStrings = {'predN'};
                varNames = myDistTable.Properties.VariableNames;
                myPos = find(ismember(varNames,myNanStrings));
                nanifyVars = varNames(myPos);
                for varCounter = 1 : length(nanifyVars)
                    myDistTable.(nanifyVars{varCounter}) = nan(nRows,1);
                end
                myNanStrings = {'predProbs'};
                varNames = myDistTable.Properties.VariableNames;
                myPos = find(ismember(varNames,myNanStrings));
                nanifyVars = varNames(myPos);
                myDistTable{1 : nRows,nanifyVars} = {nan};            
           end
           obj.distTable = myDistTable;           
           
           obj.gofMn = '';
           obj.gofSD = '';
           obj.gofBin = '';
           obj.gofDist = '';
           obj.gof = '';
      end

      function obj = getSimData(obj, myWorksheet)
           % Get all of the simulation datapoints for the VPs that will
           % be needed for calculating population statistics
           %
           % ARGUMENTS
           %  (self):      Note that the following properties should be
           %               initialized (experimental data) before calling this
           %               method:
           %                mnSDTable
           %                binTable
           %                distTable
           %  myWorksheet: A worksheet with the simulation results.
           %
           % RETURNS
           %  (self):      The VPop object is returned, but with an updated
           %               simData property.
           %
           %
           continueFlag = true;
           mnsdDataFlag = false;
           binDataFlag = false;
           distDataFlag = false;
           myMnSdData = obj.mnSDTable;
           myBinTable = obj.binTable;
           myDistTable = obj.distTable;
           % We want the {expVarID, element, elementType, time} sets
           % from the experimental summaries that will be paired
           % with real data.
           rowInfoNames = {'expVarID','interventionID','elementID','elementType','time'};
           myDataSource = cell2table(cell(0,5), 'VariableNames', rowInfoNames);
           myDataSource = '';
           if ~isempty(myMnSdData)
               myDataSource = [myDataSource; myMnSdData(:,rowInfoNames)];
               mnsdDataFlag = true;
           end
           if ~isempty(myBinTable)
               myDataSource = [myDataSource; myBinTable(:,rowInfoNames)];
               binDataFlag = true;
           end           
           if ~isempty(myDistTable)
               myDataSource = [myDataSource; myDistTable(:,rowInfoNames)];
               distDataFlag = true;
           end              
           [nrows,ncols] = size(myDataSource);
           if nrows < 1
               warning(['No mnSDTable, binTable, or distTable assigned before calling getSimData in ',mfilename,'. The data to get is not known. Exiting.'])
               continueFlag = false;
           end
           if continueFlag
               % Eliminate duplicate rows
               myDataSource = unique(myDataSource, 'rows');
               myDataSource = sortrows(myDataSource,{'interventionID','elementID','time','elementType'},{'ascend','ascend','ascend','ascend'});

               [nEntries, ~] = size(myDataSource);           

               % Info on each row of simvalues
               rowInfo = table2cell(myDataSource);
               % We will need to use a structure, like we did for results,
               % since some of the VP names will certainly be longer
               % than what is permitted by variable name lengths
               % Also, use the values
               % Row indices for the data are also added to speed subsequent calculations.

               vpIDs = getVPIDs(myWorksheet);
               nVPs = length(vpIDs);
               % Simulation results for each VP for each row
               dataValues = nan(nEntries, length(vpIDs));
               interventionIDIndex = find(ismember(rowInfoNames,'interventionID'));
               elementIDIndex = find(ismember(rowInfoNames,'elementID'));
               elementTypeIndex = find(ismember(rowInfoNames,'elementType'));
               expTimeIndex = find(ismember(rowInfoNames,'time'));
               simData.Data = dataValues;
               simData.rowInfoNames = rowInfoNames;
               simData.rowInfo = rowInfo;
               simData.vpIDs = vpIDs;
               interventionIDs = getInterventionIDs(myWorksheet);
               mnSDRows = nan(nEntries,1);
               binRows = nan(nEntries,1);
               distRows = nan(nEntries,1);

               for rowCounter = 1 : nEntries
                    interventionID = simData.rowInfo{rowCounter,interventionIDIndex};
                    elementID = simData.rowInfo{rowCounter,elementIDIndex};
                    elementType = simData.rowInfo{rowCounter,elementTypeIndex};
                    wshInterventionIndex = find(ismember(interventionIDs,interventionID));
                    expTime = simData.rowInfo{rowCounter,expTimeIndex};
                    % To avoid re-seatching for the right rows on every
                    % iteration mapel, we provide the indices here 
                    if mnsdDataFlag
                        temp = find((ismember(myMnSdData{:,'interventionID'},interventionID)) & ((myMnSdData{:,'time'})==expTime) & (ismember(myMnSdData{:,'elementID'},elementID)) & (ismember(myMnSdData{:,'elementType'},elementType)));
                        if ~isempty(temp)
                            mnSDRows(rowCounter) = temp;
                        end
                    end
                    if binDataFlag
                        temp = find((ismember(myBinTable{:,'interventionID'},interventionID)) & ((myBinTable{:,'time'})==expTime) & (ismember(myBinTable{:,'elementID'},elementID)) & (ismember(myBinTable{:,'elementType'},elementType)));
                        if ~isempty(temp)
                            binRows(rowCounter) = temp;
                        end
                    end
                    if distDataFlag
                        temp = find((ismember(myDistTable{:,'interventionID'},interventionID)) & ((myDistTable{:,'time'})==expTime) & (ismember(myDistTable{:,'elementID'},elementID)) & (ismember(myDistTable{:,'elementType'},elementType)));
                        if ~isempty(temp)
                            distRows(rowCounter) = temp;
                        end
                    end                    
                    % I don't see how to avoid looping this, given the way
                    % the data is structured.  Luckily we just get the data
                    % once
                    for vpCounter = 1 : nVPs
                        % We need to verify the desired result for the VP is
                        % present, otherwise we report the simData result as
                        % nan
                        if length(myWorksheet.results) > 0
                            curResult = myWorksheet.results{wshInterventionIndex, vpCounter};
                            % Results should be stored in a structure, we 
                            % assume if a structure is provided then it is a 
                            % valid result
                            if strcmp(class(curResult),'struct')
                                curResult = myWorksheet.results{wshInterventionIndex, vpCounter};
                                curTimeIndex = find(ismember(curResult.Names,'time'));
                                curVarIndex = find(ismember(curResult.Names,elementID));
                                curTime = curResult.Data(:,curTimeIndex);
                                curVar = curResult.Data(:,curVarIndex);
                                % Interpolate to get the time exactly right
                                interpolateValue = interp1(curTime,curVar,expTime,'linear');
                                simData.Data(rowCounter, vpCounter) = interpolateValue;
                            else
                                simData.Data(rowCounter, vpCounter) = nan;
                            end
                        else
                            simData.Data(rowCounter, vpCounter) = nan;
                        end
                    end
               end
               simData.binRows = binRows;
               simData.mnSDRows = mnSDRows;
               simData.distRows = distRows;
               obj.simData = simData;
           end
      end

      function obj = addDistTableSimVals(obj)   
          % Here we simply get sim values that are
          % fixed during optimization and add them to the 
          % dist table.  This is done at initialization to speed
          % execution.
          %
          % ARGUMENTS
          %  (self):      Note that the following properties should be
          %               initialized (experimental and simulation data) 
          %               before calling this method:
          %                distTable
          %                simData
          %
          % RETURNS
          %  (self):      The VPop object is returned, but with updated
          %               properties:
          %                distTable	
          

          myDistTable = obj.distTable;          		  
          
          mySimData = obj.simData.Data;
          
          mySimRowInfo = obj.simData.rowInfo;
          
          mySimColNames = obj.simData.rowInfoNames;
          distRowsTarget = obj.simData.distRows;
          distRowsSource = find(distRowsTarget>0); 		
          
          simInterventionIDCol = find(ismember(mySimColNames, 'interventionID'));
          simTimeCol = find(ismember(mySimColNames, 'time'));
          simElementIDCol = find(ismember(mySimColNames, 'elementID'));
          simElementTypeCol = find(ismember(mySimColNames, 'elementType'));                        
		  
          if ~isempty(distRowsSource)
              distRowsTarget = distRowsTarget(distRowsSource);
          
              % 2 step assignment to speed execution
              % first to matrix, then convert back to table.
              [nDistRows, nDistCols] = size(myDistTable);          
				
              curSimValues = nan(nDistRows,  size(mySimData,2));
              curSimValues(distRowsTarget,:) = (mySimData(distRowsSource, :));
              
              [curSimValues, I] = sort(curSimValues, 2, 'ascend');
              for rowCounter = 1 : nDistRows
                  % There should not be nan's here (without a dropout simfilter)
                  keepIndices = find(~isnan(curSimValues(rowCounter,:)));
				  
				  
                  curVals = curSimValues(rowCounter,keepIndices);                 
                  myDistTable.('predSample'){rowCounter} = curVals(keepIndices);
				  myDistTable.('predIndices'){rowCounter} = I(rowCounter,keepIndices);
                  % Also get the indices to align the exp and sim samples
                  sample1 = myDistTable.('expSample'){rowCounter};
                  [sample1Ind, sample2Ind, SC] = alignSamples(sample1, curVals);
                  myDistTable.('expCombinedIndices'){rowCounter} = sample1Ind;
                  myDistTable.('simCombinedIndices'){rowCounter} = sample2Ind;
                  myDistTable.('combinedPoints'){rowCounter} = SC;
              end
              obj.distTable = myDistTable;												
          end   

      end

      function obj = addPredTableVals(obj)
          % Once we have read the simData and generated a pw vector,
          % We need to add the predicted equivalents to the 
          % tables used for statistical comparison.
          %
          % ARGUMENTS
          %  (self):      Note that the following properties should be
          %               initialized (experimental and simulation data) 
          %               before calling this method:
          %                mnSDTable
          %                binTable
          %                distTable
          %                simData
          %
          % RETURNS
          %  (self):      The VPop object is returned, but with updated
          %               properties:
          %                mnSDTable
          %                binTable
          %                distTable
          
          myMnSdData = obj.mnSDTable;
          myBinTable = obj.binTable;
          myDistTable = obj.distTable;          
          mySimData = obj.simData.Data;
          mySimRowInfo = obj.simData.rowInfo;
          myPWs = obj.pws;
          mySimColNames = obj.simData.rowInfoNames;
                    
          mnSDRowsTarget = obj.simData.mnSDRows;
          mnSDRowsSource = find(mnSDRowsTarget>0);          
          
          binRowsTarget = obj.simData.binRows;
          binRowsSource = find(binRowsTarget>0);
                    
          distRowsTarget = obj.simData.distRows;
          distRowsSource = find(distRowsTarget>0);          
                    

          if obj.useEffN
               curN = 1/sum(myPWs.^2);
          else
              % We could use the PW cutoff here, but it seems this  
              % encourages the optimizer to try to push the weight onto a 
              % few VPs to decrease N.  Instead, let's use the number of 
              % VPs for the purpose of statistical comparison, especially
              % during optimization.
              % curN = sum(myPWs >= obj.pwCutoff);
              curN = length(myPWs);
          end

          simInterventionIDCol = find(ismember(mySimColNames, 'interventionID'));
          simTimeCol = find(ismember(mySimColNames, 'time'));
          simElementIDCol = find(ismember(mySimColNames, 'elementID'));
          simElementTypeCol = find(ismember(mySimColNames, 'elementType'));
          
          
          if ~isempty(mnSDRowsSource)
              mnSDRowsTarget = mnSDRowsTarget(mnSDRowsSource);
          
              mnSDRows = obj.simData.mnSDRows;
              mnSDRows = mnSDRows(find(mnSDRows>0));

              %  myExpMnSdData has colnames
              % {'time', 'interventionID', 'elementID', 'elementType', 'expDataID', 'expTimeVarID', 'expVarID', 'weight', 'expN', 'expMean', 'expSD', 'simN', 'simMean', 'simSD'}
              % We want to match to keep:
              [nMnSdRows, nMnSdCols] = size(myMnSdData);          

              % 2 step assignment to speed execution
              % first to matrix, then convert back to table.
              curMean = nan(nMnSdRows,1);
              curSD = nan(nMnSdRows,1);
              curSimValues = nan(nMnSdRows,length(myPWs));
              curSimValues(mnSDRowsTarget,:) = (mySimData(mnSDRowsSource, :));

              for rowCounter = 1 : nMnSdRows
                   % Note that if any simulation result values are NaN, then 
                   % these summary stats will be nan as well.  So, the default 
                   % behavior is that failed simulation runs will preclude 
                   % calculation of population statistics, which seems a 
                   % desirable failure mode.
                    curMean(rowCounter) = wtdMean(curSimValues(rowCounter,:),myPWs);
                    curSD(rowCounter) = wtdStd(curSimValues(rowCounter,:),myPWs);
              end

              myMnSdData.('predN') = (curN * ones(nMnSdRows,1));
              myMnSdData.('predMean') = (curMean);
              myMnSdData.('predSD') = (curSD);  
              obj.mnSDTable = myMnSdData;
          end
          if ~isempty(binRowsSource)
              binRowsTarget = binRowsTarget(binRowsSource);
          
              % 2 step assignment to speed execution
              % first to matrix, then convert back to table.
              [nBinRows, nBinCols] = size(myBinTable);          
              curProbs = nan(nBinRows,4);     

              curSimValues = nan(nBinRows,length(myPWs));
              curSimValues(binRowsTarget,:) = (mySimData(binRowsSource, :));
              binEdgeValues = [myBinTable{:,'binEdge1'},myBinTable{:,'binEdge2'},myBinTable{:,'binEdge3'}];          

              for rowCounter = 1 : nBinRows
                   % 'binEdge1','binEdge2','binEdge3'
                   curProbs(rowCounter,:) = wtdBinProb(curSimValues(rowCounter,:), myPWs, binEdgeValues(rowCounter,:));
              end
              myBinTable.('predN') = (curN * ones(nBinRows,1));
              myBinTable.('predBin1') = (curProbs(:,1));
              myBinTable.('predBin2') = (curProbs(:,2));
              myBinTable.('predBin3') = (curProbs(:,3));
              myBinTable.('predBin4') = (curProbs(:,4));      
              obj.binTable = myBinTable;
          end
          if ~isempty(distRowsSource)
              [nDistRows, nDistCols] = size(myDistTable); 
			  assignPWs = cell(nDistRows,1);			  
              for rowCounter = 1 : nDistRows
                  % We have already found these
                  keepIndices = myDistTable.('predIndices'){rowCounter};     
				  curPWs = myPWs(keepIndices); 														  
				  assignPWs{rowCounter} = curPWs;
              end
			  myDistTable.('predN') = (curN * ones(nDistRows,1));
			  myDistTable.('predProbs') = (assignPWs);
              obj.distTable = myDistTable;
          end      
      end

      function obj = VPop()
          % This is the constructor method for an instance of a VPop
          % (virtual population) object.
          obj.indexTable = '';
          obj.binEdges = '';
          obj.binMidPoints = '';
          obj.binProbs = [];
          obj.pws = [];
          obj.mnSDTable = '';
          obj.binTable = '';
          obj.distTable = '';          
          obj.expData = '';          
          obj.simData = '';
          obj.gofMn = '';
          obj.gofSD = '';
          obj.gofBin = '';
          obj.gofDist = '';          
          obj.gof = '';          
          obj.spreadOut = 0;
          obj.useEffN = false;
          obj.exactFlag = true; 
          obj.optimizeTimeLimit = 10*60;
          obj.optimizeType = 'pso';    
          obj.optimizePopSize = 1000;      
          obj.intSeed = -1;
          obj.tol = 1E-3;
          obj.nIters = 10000;
          obj.minEffN = 0;
      end

end

end

