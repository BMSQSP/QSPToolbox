function wsd = wtdStd(measVals,pws)
% Calculates a weighted standard deviation. Note this is a quantity
% without universal agreement on a definition from statisticians.
% So here, we adopt the definition used in NIST software. See:
% http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
%
% ARGUMENTS
%  measVals: A 1 X nVP (or nVP X 1) vector of measurement values
%  pws:      A 1 X nVP (or nVP X 1) vector of prevalence weights
%
% RETURNS
%  wsd:      weighted standard deviation
%

continueFlag = true;
wsd = nan;
% We just check if matrix/vector dimensions are OK here.
if nargin > 2
    continueFlag = false;
    warning(['A 1 x nVP (or nVP X 1) vector of values and weights are both required for ',mfilename,'. Too many input arguments.  Returning NaN.'])
elseif nargin > 1
    continueFlag = true;
else
    continueFlag = false;
    warning(['A 1 x nVP (or nVP X 1) vector of values and weights are both required for ',mfilename,'. Insufficient input arguments.  Returning NaN.'])
end

if continueFlag
    [size1r, size1c] = size(measVals);
    [size2r, size2c] = size(pws);
    if ((size1r~= size2r) || (size1c~= size2c))
        continueFlag = false;
    end
    if (min(size1r, size1c) ~= 1) || (min(size2r, size2c) ~= 1)
        continueFlag = false;
    end
    if ~continueFlag
        warning(['A 1 x nVP (or nVP X 1) vector of values and weights are both required for ',mfilename,'. Returning NaN.'])
    end
end

if continueFlag
	% We allow an error of eps * the number of VPs * 10
	% We could in strictly enforce eps * the number of VPs but
	% this criterion is relaxed a little.
	% 10^-14 was used as a minimum in the bin probabilities.  The
	% bin probabilities should all be fixed and
	% PWs should all be corrected to sum to 1, but
	% here 10^-14 is set as a minimum allowed error in the check
	numericTolerance = max(eps(1) * length(pws) * 10,1E-14);	
    if ((sum(pws) > (1 + numericTolerance)) || (sum(pws) < (1 - numericTolerance)))
        warning(['Supplied pws to ',mfilename,' should sum to 1. Returning NaN.'])
        continueFlag = false;
    end
end      
    
if continueFlag
    wMean = wtdMean(measVals,pws);
    wVar = sum(pws.*((measVals - wMean).^2));
    % Note that for n in the wsd, we use the number of pws above a numeric
    % threshold (here, 1 E -16 as was done in the MAPEL paper).   
    % If an effective N is used for statistical comparisons, the effective
    % N is implemented at the time of the statistical testing (e.g. F-test,
    % t-test, ...)
    n = sum(pws > numericTolerance);
    % Use the unbiased correction so it gives the
    % same answer as R's SD function
    wVar = wVar*n/(n-1);
    wsd = sqrt(wVar);  
    % infinity can be returned
    % with a PW of ~1 coupled with
    % very small PWs.  The std
    % of 1 value should be 0.
    if isinf(wsd)
        wsd = 0;
    end
end
end