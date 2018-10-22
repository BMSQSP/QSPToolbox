function copiedWorksheet = copyWorksheet(myWorksheet, vpIDs, flagCopyResults, flagCheckVPIDs)
% Make a copy of a worksheet; if VPs are specified, they will be selected.
%
% ARGUMENTS
% myWorksheet:     A worksheet
% vpIDs:           A cell array, 1XnVP, of string identifiers
% flagCopyResults: Whether to copy the results from the current worksheet.
%                   Included as a potential memory saver.
% flagCheckVPIDs:  Whether to check the VPIDs.  Default is true, but
%                   requires more time especially with large worksheets
%
% RETURNS
% copiedWorksheet
%

% Perform initial checks on the provided arguments
flagContinue = false;
copiedWorksheet = createWorksheet();

if nargin > 4
    warning([mfilename,' requires input argument: myWorksheet, and optionally vpIDs, flagCopyResults, flagCheckVPIDs.  Too many arguments provided.'])
elseif nargin == 0 
    warning([mfilename,' requires input argument: myWorksheet, and optionally vpIDs, flagCopyResults, flagCheckVPIDs.  Insufficient arguments provided.'])
elseif nargin == 1 
    vpIDs = getVPIDs(myWorksheet);
    flagCopyResults = true;
	% If we're not given VPIDs we won't check them!
	flagCheckVPIDs = false;
    flagContinue = true;
elseif nargin == 2 
    flagCopyResults = true;
	flagCheckVPIDs = true;
    flagContinue = true;
elseif nargin == 3 
    flagContinue = true; 
	flagCheckVPIDs = true;	
elseif nargin == 4 
    flagContinue = true;     
end

if flagContinue
    vpIDsWorksheet = getVPIDs(myWorksheet);
	if flagCheckVPIDs
		% The provided vpIDs may not contain the same members as the worksheet.
		remainingOriginalIndices = (1 : length(vpIDsWorksheet))';
		vpOriginalWorksheetIndices = nan(1,length(vpIDs));
		for vpCounter = 1 : length(vpIDs)
			vpID = vpIDs{vpCounter};
			oldIndex=find(ismember(vpIDsWorksheet, vpID));
			if length(oldIndex) > 1
				warning(['Degenerate VP IDs in worksheet provided to ',mfilename,'.'])
				flagContinue = false;
			elseif length(oldIndex) < 1
				warning(['Missing VP IDs in worksheet provided to ',mfilename,'.'])
				flagContinue = false;
			else
				vpOriginalWorksheetIndices(1,vpCounter)=remainingOriginalIndices(oldIndex(1));
				remainingOriginalIndices(oldIndex(1)) = [];
				vpIDsWorksheet(oldIndex(1))=[];
			end
		end
	else
		vpOriginalWorksheetIndices = find(ismember(vpIDsWorksheet,vpIDs));
	end	
end
    
% Copy the worksheet fields over
if flagContinue    
    fieldNames = fields(myWorksheet);    
    for fieldCounter = 1 : length(fieldNames)
        curField = fieldNames{fieldCounter};
        if strcmp(curField,'vpDef')
            copiedWorksheet.(curField) = myWorksheet.(curField)(vpOriginalWorksheetIndices);
        elseif strcmp(curField,'axisProps')
            copiedWorksheet.axisProps = struct();
            copiedWorksheet.(curField).axisDef = myWorksheet.(curField).axisDef;
            myCoefficients = getVPCoeffs(myWorksheet);
            copiedWorksheet.(curField).axisVP = myWorksheet.(curField).axisVP;
            copiedWorksheet.(curField).axisVP.coefficients = myCoefficients(:,vpOriginalWorksheetIndices);
        elseif strcmp(curField,'responseTypes')
            % Since RTEs are objects, the copy is written down to the
            % level of RTEs.  This isn't strictly necessary,
            % but may help in the future if we want to use
            % memory references for the objects.
            [nRTs, dummy] = size(myWorksheet.(curField));
            copiedWorksheet.(curField) = {};
            for rtCounter = 1 : nRTs
                copiedWorksheet.(curField){rtCounter,1} = struct();
                rtFields = fields(myWorksheet.(curField){rtCounter,1});
                for rtFieldCounter = 1 : length(rtFields)
                    rtField = rtFields{rtFieldCounter};
                    if ~(strcmp(rtField,'elements')) 
                        copiedWorksheet.(curField){rtCounter,1}.(rtField) = myWorksheet.(curField){rtCounter,1}.(rtField);
                    else
                        curElements = myWorksheet.(curField){rtCounter,1}.(rtField);
                        [nElements, dummy] = size(curElements);
                        myWorksheet.(curField){rtCounter,1}.(rtField) = {};
                        for eCounter = 1 : nElements
                            copiedWorksheet.(curField){rtCounter,1}.(rtField){eCounter,1} = curElements{eCounter,1};
                        end
                    end
                end
            end            
        elseif strcmp(curField,'results')
            if ~(flagCopyResults)
                copiedWorksheet.(curField) = cell(0,0);
            else
                [nInterventionResults, nVPresults] = size(myWorksheet.(curField));
                if (nVPresults > 0)
                    copiedWorksheet.(curField) = myWorksheet.(curField)(:,vpOriginalWorksheetIndices);
                else
                    copiedWorksheet.(curField) = {};
                end
            end
        else
            % We directly "copy" here with "=": interventions, expData, simProps,
            % variantProps, model, compiled
            copiedWorksheet.(curField) = myWorksheet.(curField);
        end
    end
else
    warning(['Could not complete ',mfilename,'. Returning original worksheet.'])
    copiedWorksheet = myWorksheet;
end

end