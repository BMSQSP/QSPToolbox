classdef linearCalibrateOptions
% Here, we define the linearCalibrateOptions class to run a VPopOptimization
%
%  Properties
%   verbose: 				(boolean, optional): Whether to print real-time status to 
% 				    		command window. Default is 'false'
%   cdfProbsToFit: 			(vector, optional but recommended): Probabilities of the cumulative
% 							distribution functions (CDF) to fit. Default is to fit all points in the
% 							CDF, which may take a long time.
%  toleranceScalingFactor:  (scalar, optional): Scaling factor for the optimization tolerance 
%							(e.g., '0.1' would decrease the tolerance by 10).  Default is 1.
%  binTableGroupWeight: 				(scalar) Weight for data in myVPop.binTable
%  distTableGroupWeight: 				(scalar) Weight for data in myVPop.distTable
%  brTableRECISTGroupWeight:			(scalar) Weight for data in myVPop.brTableRECIST
%  rTableRECISTGroupWeight: 			(scalar) Weight for data in myVPop.rTableRECIST
%  mnSDTableGroupWeight:				(scalar) Weight for data in myVPop.mnSDTable
%										GroupWeights Specify fitting weights for each of the
% 										different types of data. Default is to set all of the weights equal to 1,
% 										so that all of the different data groups are weighed equally.
%  priorPrevalenceWeightAssumption:		(string, optional) A priori
% 										knowledge of prevalence weights is required to renormalize VP prevalence weights
% 										in cases where some VP simulation values are NaN since they dropped off
% 										therapy. With VPs missing, the prevalence weights need to be renormalized
% 										before calculating the weighted mean and standard deviation of a
% 										simulation variable across the virtual population. Three options are
% 										available: 
%										(1) 'uniform': assume uniform prevalence weights (default) assumes that
% 										all VPs have the same prevalence weight -- this isn't expected to result
% 										in inaccurate fitting if the probability of dropping off therapy doesn't
% 										correlate with prevalence weight, since then the average prevalence
% 										weight of a large chunk of NaN VPs should be approximately the uniform
% 										prevalence weight, 1/nVP, and for a small chunk of NaN VPs, the
% 										prevalence weight renormalization correction factor should be negligible;
% 										(2) 'specified': assume prevalence weights uses the prevalence weights
% 										specified in myVPop.pws as the prior prevalence weights; and 
%										(3)  'ignoreDropout': do not fit simulation data containing NaN VPs, ignores any data in which
% 										at least some of the VPs are NaN which can happen with modeled dropouts
%  responseValTransformation:	(string, optional): Specifies whether the
% 								observations should be transformed  so that the residuals relative to the
% 								values of the observations are what are being fitting. Options:
% 								'relative' (default) and 'none'
% 	createPlotsOfOptimizedVPop: (boolean) whether to generate plots of the fit
% 	createPlotsOfInputVPop:		(boolean) whether to generate plots of the
% 								input myVPop. If 'true', this assumes that myVPop is a calibrated VPop.
%
   properties
      verbose
      cdfProbsToFit
	  toleranceScalingFactor
	  binTableGroupWeight
	  distTableGroupWeight
	  brTableRECISTGroupWeight
	  rTableRECISTGroupWeight
      mnSDTableGroupWeight
	  priorPrevalenceWeightAssumption
	  responseValTransformation
	  createPlotsOfOptimizedVPop
	  createPlotsOfInputVPop
	  
   end
   methods
      function obj = set.verbose(obj,myValue)
          if (myValue == true) || (myValue == false)
              obj.verbose = myValue;
          else
            error(['Invalid verbose value in ',mfilename,', expecting true/false'])
          end
      end
      function obj = set.cdfProbsToFit(obj,myValue)
          if isnumeric(myValue) 
              obj.cdfProbsToFit = myValue;
          else
            error(['Invalid cdfProbsToFit value in ',mfilename,', expecting a numeric vector'])
          end
      end 	  
      function obj = set.toleranceScalingFactor(obj,myValue)
          if isnumeric(myValue) 
              obj.toleranceScalingFactor = myValue;
          else
            error(['Invalid toleranceScalingFactor value in ',mfilename,', expecting a numeric value'])
          end
      end 	  	  
      function obj = set.binTableGroupWeight(obj,myValue)
          if isnumeric(myValue) 
              obj.binTableGroupWeight = myValue;
          else
            error(['Invalid binTableGroupWeight value in ',mfilename,', expecting a numeric value'])
          end
      end	
      function obj = set.distTableGroupWeight(obj,myValue)
          if isnumeric(myValue) 
              obj.distTableGroupWeight = myValue;
          else
            error(['Invalid distTableGroupWeight value in ',mfilename,', expecting a numeric value'])
          end
      end	  	  
      function obj = set.brTableRECISTGroupWeight(obj,myValue)
          if isnumeric(myValue) 
              obj.brTableRECISTGroupWeight = myValue;
          else
            error(['Invalid brTableRECISTGroupWeight value in ',mfilename,', expecting a numeric value'])
          end
      end	
      function obj = set.rTableRECISTGroupWeight(obj,myValue)
          if isnumeric(myValue) 
              obj.rTableRECISTGroupWeight = myValue;
          else
            error(['Invalid rTableRECISTGroupWeight value in ',mfilename,', expecting a numeric value'])
          end
      end	
      function obj = set.priorPrevalenceWeightAssumption(obj,myValue)
          if (strcmpi(myValue,'uniform') ||...
               strcmpi(myValue,'specified') ||...
               strcmpi(myValue,'ignoreDropout'))
            obj.priorPrevalenceWeightAssumption = myValue;
          else
            error(['Invalid priorPrevalenceWeightAssumption setting in ',mfilename,', allowed setting: "uniform," "specified," "ignoreDropout."'])
          end
      end
      function obj = set.responseValTransformation(obj,myValue)
          if (strcmpi(myValue,'relative') ||...
               strcmpi(myValue,'none'))
            obj.responseValTransformation = myValue;
          else
            error(['Invalid responseValTransformation setting in ',mfilename,', allowed setting: "relative," "none."'])
          end
      end
      function obj = set.createPlotsOfOptimizedVPop(obj,myValue)
          if (myValue == true) || (myValue == false)
              obj.createPlotsOfOptimizedVPop = myValue;
          else
            error(['Invalid createPlotsOfOptimizedVPop value in ',mfilename,', expecting true/false'])
          end
      end	  
      function obj = set.createPlotsOfInputVPop(obj,myValue)
          if (myValue == true) || (myValue == false)
              obj.createPlotsOfInputVPop = myValue;
          else
            error(['Invalid createPlotsOfInputVPop value in ',mfilename,', expecting true/false'])
          end
      end	  

          
      
      function value = get(obj,propName)
          switch propName
              case 'verbose'
                  value = obj.verbose;
              case 'cdfProbsToFit'
                  value = obj.cdfProbsToFit;				  
              case 'toleranceScalingFactor'
                  value = obj.toleranceScalingFactor;                  
              case 'binTableGroupWeight'
                  value = obj.binTableGroupWeight;
              case 'distTableGroupWeight'
                  value = obj.distTableGroupWeight;                  
              case 'brTableRECISTGroupWeight'
                  value = obj.brTableRECISTGroupWeight;
              case 'rTableRECISTGroupWeight'
                  value = obj.rTableRECISTGroupWeight;
              case 'mnSDTableGroupWeight'
                  value = obj.mnSDTableGroupWeight;
              case 'priorPrevalenceWeightAssumption'
                  value = obj.priorPrevalenceWeightAssumption;
              case 'responseValTransformation'
                  value = obj.responseValTransformation;
              case 'createPlotsOfOptimizedVPop'
                  value = obj.createPlotsOfOptimizedVPop;
              case 'createPlotsOfInputVPop'
                  value = obj.createPlotsOfInputVPop; 				  
              otherwise
                  error(['Error: ',propName ,' is not a valid ',mfilename,' property.'])
          end
      end      

      % The constructor method must have the same name as the class
      function obj = linearCalibrateOptions()
			obj.verbose = false;
		  
			obj.cdfProbsToFit = 0.05:0.05:0.95;

			% Set as less than 1 to reduce tolerance for optimization
			obj.toleranceScalingFactor = 1;

			% Weights of different types of data
			obj.binTableGroupWeight = 1;
			obj.distTableGroupWeight = 1;
			obj.brTableRECISTGroupWeight = 1;
			obj.rTableRECISTGroupWeight = 1;
			obj.mnSDTableGroupWeight = 1;


			obj.priorPrevalenceWeightAssumption = 'uniform';

			% Options:
			% 'relative' (default)
			% 'none'
			obj.responseValTransformation = 'relative';

			obj.createPlotsOfInputVPop = false;
			obj.createPlotsOfOptimizedVPop = false;		  
			
      end
   end
end