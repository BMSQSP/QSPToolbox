function myExpDataTable = convertMapelOptionsToExpDataTable(myWorksheet, myResponseTypeID)
% This function takes experimental data embedded in a worksheet and
% converts it to a data table format for use in MAPEL, guided by a
% response type
%
% ARGUMENTS:
% myWorksheet
% myResponseTypeID
%
% RETURNS
% myExpDataTable

continueFlag = true;
if nargin > 2
    continueFlag = false;
    warning(['Too many input arguments for ',mfilename,'. Should provide: myWorksheet, myResponseTypeID.'])
    continueFlag = false;
elseif nargin > 1
    continueFlag = true;
else
    warning(['Insufficient input arguments for ',mfilename,'. Should provide: myWorksheet, myResponseTypeID.'])
    continueFlag = false;
end

if continueFlag
    allResponseTypeIDs = getResponseTypeIDs(myWorksheet);
    
    if sum(ismember(allResponseTypeIDs, myResponseTypeID)) < 1
        warning(['Specified responseTypeID not found in the provided worksheet in call to ',mfilename,'.'])
        continueFlag = false;
    
    else
        myResponseType = getResponseType(myWorksheet, myResponseTypeID);
        [nRTE,~] = size(myResponseType.elements);
        for responseTypeElementCounter = 1 : nRTE
            myRTE = myResponseType.elements{responseTypeElementCounter};
            if ~(strcmp(class(myRTE), 'responseTypeElementPoints'))
                warning(['Included responseTypeElementID ,', myRTE.id,', not of class responseTypeElementPoints, this case is not currently handled in ',myResponseTypeID,'.'])
                continueFlag = false;                
            end
        end
    end
end

if continueFlag
    [nRTE, dummy] = size(myResponseType.elements);
    for responseTypeElementCounter = 1 : nRTE
        % We may want to add a filter for the first data point,
        % in case initial conditions are specified,
        % or this can be done post-processing.
        myRTE = myResponseType.elements{responseTypeElementCounter};
        myInterventionID = myRTE.interventionID;
        myElementID = myRTE.modelYVar;
        myElementType = myRTE.modelYVarType;
        myExpDataID = myRTE.expDataID;
        myExpTimeVarID = myRTE.expDataTimeVar;
        myExpVarID = myRTE.expDataYVar;
        if responseTypeElementCounter == 1
            myExpDataTable = createExpDataTable(myWorksheet, myInterventionID, myElementID, myElementType, myExpDataID, myExpTimeVarID, myExpVarID);
        else
            myExpDataTable = createExpDataTable(myWorksheet, myInterventionID, myElementID, myElementType, myExpDataID, myExpTimeVarID, myExpVarID, myExpDataTable);
        end
    end
else
    warning(['Unable to complete ',mfilename,', exiting.'])
end    