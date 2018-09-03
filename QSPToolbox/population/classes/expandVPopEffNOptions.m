classdef expandVPopEffNOptions
% Options and data required to to run expandVPopEffN
%
%  suffix:             (Required) a text descriptor string that will be included in what 
%                      will be written to file.  This is also used
%                      in setting VP identities.
%                       Default = '' and should be specified  
%  wsIterCounter:      (Optional) tracks the iterations through the algorithm.  WARNING: Keep
%                      incrementing to avoid issues with repeated VPIDs, and do not repeat
%                      a previous value.
%                       Default = 0  
%  targetEffN:         (Optional) algorithm stops when it is achieves a valid VPop
%                      with this effN
%                       Default = 100  
%  maxNewPerIter:      (Optional) Maximum new VPs we can add per iteration.  Set to 
%                      -1 to use the VPop effN
%                       Default = 50  
%  expandCohortSize:   (Optional) size of the cohort to generate for testing when
%                      making new VPs
%                       Default = 2500                       
%  effNDelta:          (Optional) How much to increment the minEffN each iteration
%                       Default = 2
%  hotStartVPop:       (Optional) If we already have a good VPop for mapelOptions and either the 
%                      current or a previous worksheet iteration, use this here.  
%					   It will just be important that the first VPs in the worksheet
%                      haven't been rearranged since the VPop was found.
%                      That is, the worksheet may have new/more VPs, but the first N VPs
%                      in the worksheet corresponse to the N VPs in the VPop used in the fit. 
%                      It will be used in the initial worksheet iteration in case the 
%                      algorithm has problems.  Enter '' if there is none.
%                       Default = ''
%  minPVal:            (Optional) minimum p-value for a good VPop
%                       Default = 0.9
%  nTries:             (Optional) Number of times to try a solution before triggering expansion.
%                       Default = 5
%  nRetries:           (Optional) Number of times to try a solution that is not quite at acceptance.
%                       Default = 30
%  restartPVal:        (Optional) p-value that will trigger a restart of the algorithm.
%                       Default = 1E-4
%  expandRandomStart:  (Optional) If 0, the prevalence 
%                      weights are not perturbed from initial guesses on the first iteration or
%                      following a successful solution for the population and expansion.
%                      Otherwise the transformed initial bin
%                      probabilities are perturbed 
%                      with random (normally distributed) noise that has
%                      normalized standard deviation of the specified value.  
%                       Default = 0; 
%  verbose:            (Optional) report progress to screen
%                       Default = true;
% 

   properties
		suffix
		wsIterCounter
		targetEffN
		maxNewPerIter
		expandCohortSize
		effNDelta
		minPVal
		nTries
		nRetries
		restartPVal
		expandRandomStart
		verbose
   end

   methods     

      function obj = set.suffix(obj,mySuffix)
          if ischar(mySuffix)
              obj.suffix = (mySuffix);
          else
              error(['Invalid suffix specified in ',mfilename,', a string should be specified.'])
          end
      end         
   
      function obj = set.wsIterCounter(obj,myWsIterCounter)
          failFlag = false;
          if isnumeric(myWsIterCounter) 
              if isequal(size(myWsIterCounter),[1 1])
                  if round(myWsIterCounter) >= 0
                      obj.wsIterCounter = round(myWsIterCounter);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid wsIterCounter specified for ',mfilename,'. A nonnegative number should be specified.'])
          end
      end   
   
      function obj = set.targetEffN(obj,myTargetEffN)
          failFlag = false;
          if isnumeric(myTargetEffN) 
              if isequal(size(myTargetEffN),[1 1])
			  	  if mod(myTargetEffN,1) ~= 0
					  warning(['A whole number is required for targetEffN in ',mfilename,'. The provided value will be rounded.']) 
				  end
                  if round(myTargetEffN) >= 0
                      obj.targetEffN = round(myTargetEffN);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid targetEffN specified for ',mfilename,'. A whole number should be specified.'])
          end
      end      
   
      function obj = set.maxNewPerIter(obj,myMaxNewPerIter)
          failFlag = false;
          if isnumeric(myMaxNewPerIter) 
              if isequal(size(myMaxNewPerIter),[1 1])
                  if round(myMaxNewPerIter) > 0
                      obj.maxNewPerIter = round(myMaxNewPerIter);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid maxNewPerIter specified for ',mfilename,'. A positive number should be specified.'])
          end
      end     
 
      function obj = set.expandCohortSize(obj,myExpandCohortSize)
          failFlag = false;
          if isnumeric(myExpandCohortSize) 
              if isequal(size(myExpandCohortSize),[1 1])
                  if round(myExpandCohortSize) >= 0
                      obj.expandCohortSize = round(myExpandCohortSize);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid expandCohortSize specified for ',mfilename,'. A nonnegative number should be specified.'])
          end
      end    

      function obj = set.effNDelta(obj,myEffNDelta)
          failFlag = false;
          if isnumeric(myEffNDelta) 
              if isequal(size(myEffNDelta),[1 1])
				  if mod(myEffNDelta,1)~=0
					  warning(['A whole number is required for effNDelta in ',mfilename,'. The provided value will be rounded.']) 
				  end
                  if round(myEffNDelta) >= 0.5
                      obj.effNDelta = round(myEffNDelta);
                  elseif myEffNDelta >= 0
                      obj.effNDelta = 0;
					  warning(['A zero effNDelta is specified in ',mfilename,'. Note the calling function may never terminate if an increased effN is required.']) 
                  else
                      failFlag = true;					  
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid effNDelta specified for ',mfilename,'. A positive number should be specified.'])
          end
      end  	  
		
      function obj = set.minPVal(obj,myMinPVal)
          failFlag = false;
          if isnumeric(myMinPVal) 
              if isequal(size(myMinPVal),[1 1])
                  if ((myMinPVal >= 0) && (myMinPVal <= 1))
                      obj.minPVal = myMinPVal;
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid minPVal specified for ',mfilename,'. A number between 0 and 1 should be specified.'])
          end
      end  
	
      function obj = set.nTries(obj,myNTries)
          failFlag = false;
          if isnumeric(myNTries) 
              if isequal(size(myNTries),[1 1])
                  if round(myNTries) > 0
                      obj.nTries = round(myNTries);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid nTries specified for ',mfilename,'. A positive number should be specified.'])
          end
      end  	  	

      function obj = set.nRetries(obj,myNRetries)
          failFlag = false;
          if isnumeric(myNRetries) 
              if isequal(size(myNRetries),[1 1])
                  if round(myNRetries) > 0
                      obj.nRetries = round(myNRetries);
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid nRetries specified for ',mfilename,'. A positive number should be specified.'])
          end
      end  		  
	  
      function obj = set.restartPVal(obj,myRestartPVal)
          failFlag = false;
          if isnumeric(myRestartPVal) 
              if isequal(size(myRestartPVal),[1 1])
                  if ((myRestartPVal >= 0) && (myRestartPVal <= 1))
                      obj.restartPVal = myRestartPVal;
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid restartPVal specified for ',mfilename,'. A number between 0 and 1 should be specified.'])
          end
      end  	  
	  
      function obj = set.expandRandomStart(obj,myExpandRandomStart)
          failFlag = false;
          if isnumeric(myExpandRandomStart) 
              if isequal(size(myExpandRandomStart),[1 1])
                  if (myExpandRandomStart >= 0)
                      obj.expandRandomStart = myExpandRandomStart;
                  else
                      failFlag = true;
                  end
              else
                  failFlag = true;
              end
          else
              failFlag = true;
          end
          if failFlag
              error(['Invalid expandRandomStart specified for ',mfilename,'. A number >= 0 should be specified.'])
          end
      end  	     	  
	  
	  function obj = set.verbose(obj,myVerbose)
          if islogical(myVerbose)
               obj.verbose = myVerbose;
          else
               error(['Property verbose in ',milename,' must be true or false.'])
          end
      end         
 
      
      function value = get(obj,propName)
          switch propName
              case 'suffix'
                  value = obj.suffix;
              case 'wsIterCounter'
                  value = obj.wsIterCounter;
              case 'targetEffN'
                  value = obj.targetEffN;  
              case 'maxNewPerIter'
                  value = obj.maxNewPerIter;                   
              case 'expandCohortSize'
                  value = obj.expandCohortSize;				  
              case 'effNDelta'
                  value = obj.effNDelta;
              case 'minPVal'
                  value = obj.minPVal;
              case 'nTries'
                  value = obj.nTries;
              case 'nRetries'
                  value = obj.nRetries;
              case 'restartPVal'
                  value = obj.restartPVal;	
              case 'expandRandomStart'
                  value = obj.expandRandomStart;
              case 'verbose'
                  value = obj.verbose;				  
              otherwise
                  error(['Error: ',propName ,' is not a valid ',mfilename,' property.'])
          end
      end

      % TODO: ADD METHODS TO VERIFY AGAINST WORKSHEET
            
      function obj = expandVPopEffNOptions()
          % This is the constructor method for instances of 
          % expandVPopEffNOptions objects.
          obj.suffix = '';
          obj.wsIterCounter = 0;
          obj.targetEffN = 100;
          obj.maxNewPerIter = 50;          
		  obj.expandCohortSize = 2500;
          obj.effNDelta = 2;
          obj.minPVal = 0.9;
          obj.nTries = 5;
          obj.nRetries = 30;
		  obj.restartPVal = 1E-4;
		  obj.expandRandomStart = 0;
		  obj.verbose = true;		  
      end
   end
end