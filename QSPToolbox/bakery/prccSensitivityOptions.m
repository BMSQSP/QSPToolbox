classdef prccSensitivityOptions
% Options object for running PRCC sensitivity analysis.
% Note that the sensitivity analysis must be run with the returned
% worksheet from a sampling step.
%
% PROPERTIES
% analyzeElementResultIDs: Worksheet interventions to run analysis
%                           on. 
% analyzeInterventionIDs:  Worksheet interventions to run analysis
%                           on.
% analyzeTimes:            (Required) Simulation time points to run analysis
%                           at.  If not specified, the analysis will be run
%                           at all sampled timepoints.  Also,
%                           this should correstpond to an exact sampled
%                           simulation timepoint as interpolation is not
%                           implemented
% poolClose:                Whether to close the pool at the end of
%                            calculations
% poolRestart:              Whether to restart the pool at the beginning 
%                            of calculations
% clusterID:                Which cluster to use, parallel.defaultClusterProfile is default
% nWorkers:                 How many workers in the pool to use.  Set to nan to
%                            use the default (default).
%
% TODO: add in subsampling to assess if enough VPs are run
%        for consistent results
    
   properties
      analyzeElementResultIDs
      analyzeInterventionIDs
      analyzeTimes
      poolClose
      poolRestart
      clusterID
      nWorkers      
   end
   
   methods

      function obj = set.analyzeElementResultIDs(obj,myAnalyzeElementResultIDs)
          if iscell(myAnalyzeElementResultIDs)
              obj.analyzeElementResultIDs = (myAnalyzeElementResultIDs);
          else
              error(['Invalid analyzeElementResultIDs specified for ',mfilename,', a cell array of variable ID strings should be specified.'])
          end
      end            
       
      function obj = set.analyzeInterventionIDs(obj,myAnalyzeInterventionIDs)
          if iscell(myAnalyzeInterventionIDs)
              obj.analyzeInterventionIDs = (myAnalyzeInterventionIDs);
          else
              error(['Invalid analyzeInterventionIDs specified for ',mfilename,', a cell array of intervention ID strings should be specified.'])
          end
      end         
      
      function obj = set.analyzeTimes(obj,myAnalyzeTime)
          failFlag = false;
          if isnumeric(myAnalyzeTime)
              obj.analyzeTimes = (myAnalyzeTime);
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid analyzeTimes specified for ',mfilename,'.'])
          end
      end      
      
      function obj = set.poolClose(obj,myInput)
          if islogical(myInput)
               obj.poolClose = myInput;
          else
               error(['Property poolClose in ',mfilename,' must be logical.'])
          end
      end  	  
	  
      function obj = set.poolRestart(obj,myInput)
          if islogical(myInput)
               obj.poolRestart = myInput;
          else
               error(['Property poolRestart in ',mfilename,' must be logical.'])
          end
      end       
      
      function obj = set.nWorkers(obj,myNWorkers)
          if (isnumeric(myNWorkers) == true)
              obj.nWorkers = myNWorkers;
          else
            error(['Invalid nWorkers value in ',mfilename,'.'])
          end
      end      
      
      function obj = set.clusterID(obj,myClusterID)
          if ischar(myClusterID)
              obj.clusterID = myClusterID;
          else
            error(['Invalid clusterID in ',mfilename,'.'])
          end
      end          

      function value = get(obj,propName)
          switch propName                 
              case 'analyzeElementResultIDs'
                  value = obj.analyzeElementResultIDs;   
              case 'analyzeInterventionIDs'
                  value = obj.analyzeInterventionIDs;  
              case 'analyzeTimes'
                  value = obj.analyzeTimes;   
              case 'poolRestart'
                  value = obj.poolRestart;
              case 'poolClose'
                  value = obj.poolClose;
              case 'nWorkers'
                  value = obj.nWorkers;
              case 'clusterID'
                  value = obj.clusterID;                   
          end
      end     
      
      function passCheck = verify(obj,myWorksheet)
          passCheck = true;
          myOutputVariables = obj.analyzeElementResultIDs;
          mySampleTimes = obj.analyzeTimes;
          myInterventionIDs = obj.analyzeInterventionIDs;

          if ~ischar(myOutputVariables) && ~iscellstr(myOutputVariables)
              warning(['OutputVariable is not a single string or a list of strings in ',mfilename,'.  Exiting.'])
              passCheck = false;
          end
          
          allElementResultIDs = myWorksheet.simProps.saveElementResultIDs;
          if sum(~ismember(myOutputVariables,allElementResultIDs)) > 0
              warning(['myWorkSheet is missing the desired OutputVariables in ',mfilename,'.  Exiting.'])
              continueFlag = false;
          end
          
          if sum(mySampleTimes >= 0) < length(mySampleTimes)
              warning(['mySampleTimes must be a vector of nonnegative values in ',mfilename,'.  Exiting.'])
              passCheck = false;
          end
          
          allSampleTimes = myWorksheet.simProps.sampleTimes;
          if sum(ismember(mySampleTimes,allSampleTimes)) < length(mySampleTimes)
              warning(['mySampleTimes is not a subset of SampleTimes in myWorkSheet in ',mfilename,'.  Exiting.'])
              passCheck = false;
          end
          
          allInterventionIDs = getInterventionIDs(myWorksheet);
          if sum(ismember(myInterventionIDs,allInterventionIDs)) < length(myInterventionIDs)
              warning(['myInterventionIDs is not a subset of the intervention IDs in myWorkSheet in ',mfilename,'.  Exiting.'])
              passCheck = false;
          end          
          
          allInterventionIDs = getInterventionIDs(myWorksheet);
          allVPIDs = getVPIDs(myWorksheet);
          [dummy, nVPResults]= size(myWorksheet.results);
          if nVPResults ~= length(allVPIDs)
              warning(['Results for all VPs not available from worksheet in ',mfilename,'.  Exiting.'])
              passCheck = false;
          end
          
%           if (obj.nBootstraps < 1) && ~(strcomp(subSampleSplitType,'none'))
%               warning(['If splitting into subsamples in ',mfilename,', nBootstraps must be > 0.'])
%               passCheck = false;
%           end    
          
          if passCheck
              for interventionCounter = 1 : length(myInterventionIDs)
                  myInterventionID = myInterventionIDs{interventionCounter};
                  
                  interventionIndex = find(ismember(allInterventionIDs, myInterventionID));
                  checkResults = myWorksheet.results(interventionIndex,:);
                  if sum(sum(arrayfun(@(i) isstruct(checkResults{i}),1:nVPResults)) ) < length(allVPIDs)
                      %failCheckVPs = ~(arrayfun(@(i) isstruct(checkResults{i}),1:nVPResults));
                      %failCheckVPIndices = find(failCheckVPs);
                      warning(['Results for all VPs not available from worksheet in ',mfilename,'.'])
                      passCheck = false;
                  else
                      for vpCounter = 1 : nVPResults
                          if sum(ismember(checkResults{1, vpCounter}.Names,obj.analyzeElementResultIDs)) < length(obj.analyzeElementResultIDs)
                              passCheck = false;
                          end
                      end
                      if ~passCheck
                          warning(['Not all analyzeElementResultIDs identified in the worksheet results provided to ',mfilename,'.'])
                      end
                  end
              end
          end
          if passCheck
              ntimefail = 0;
              for vpCounter = 1 : nVPResults
                  timeIndex = find(ismember(checkResults{1, vpCounter}.Names,'time'));
                  timeVals = checkResults{1, vpCounter}.Data(:,timeIndex);
                  % Might want to add a tolerance here... but
                  % MATLAB has been generally robust for small
                  % numerical issues so far
                  if sum(timeVals == obj.analyzeTimes) < 1
                      ntimefail = ntimefail+1;
                  end
              end
              if ntimefail > 0
                  passCheck = false;
                  warning(['Not all sample times align with the desired analyzeTime in the results provided to ',mfilename,'.'])
                  passCheck = false;
              end
          end
      end
      
      % The constructor method must have the same name as the class
      function obj = prccSensitivityOptions() 
          obj.analyzeElementResultIDs={};
          obj.analyzeInterventionIDs={};          
          obj.analyzeTimes=0;  
          obj.poolClose = true;  
		  obj.poolRestart = true;		  
          obj.clusterID = parallel.defaultClusterProfile;
          obj.nWorkers = nan;          
      end
   end
end