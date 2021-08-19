function lmResults = regressPD2ScaleFactor(myMapelOptions)
% This function takes a mapelOptions or VPop object and returns 
% coefficients for a predicted PD2 scaling rate.
%
% ARGUMENTS
%  myMapelOptions:         (Required) a mapelOptionsRECIST or VPopRECIST
%                                     instance.
%
% Returns
%  lmResults:               Output of the regression.  This gives the
%                            a regression model for the observed PD normalized 
%                            to the number of patients recorded as nonPD2 
%                            (expN) at the first lesion scan.  Multiply by 
%                            expN to get the predicted expNPD21LS.  To
%                            given an expNPD21LS, one could use the
%                            predSD, predPD with the regression coefficients
%                            and scale by the the predN.
%                            
%


continueFlag = true;
if nargin > 1
    warning(['Too many input arguments provided to ',mfilename,'.  Requires: mapelOptionsRECIST or VPopRECIST object.'])
    continueFlag = false;
elseif nargin < 1
    continueFlag = false;
    warning(['Insufficient input arguments provided to ',mfilename,'.  Requires: mapelOptionsRECIST or VPopRECIST object.'])
end

if continueFlag
    if ~ismember(class(myMapelOptions),{'mapelOptionsRECIST','VPopRECIST'})
        continueFlag = false;
        warning(['Input not recognized in call to ',mfilename,'.  Requires: mapelOptionsRECIST or VPopRECIST object.'])
    end       
end


if continueFlag
    % First find the blocks by the PD2 measures
    testData = myMapelOptions.brTableRECIST(:,{'interventionID'});
    [C,IA,IC] = unique(testData,'rows','stable');
    testRows = find(~isnan(myMapelOptions.brTableRECIST{:,{'expNPD21LS'}}));
    regRows = intersect(testRows, IA);
    pred = table2array(myMapelOptions.brTableRECIST(regRows,{'expSD','expPD'}));
    resp = table2array(myMapelOptions.brTableRECIST(regRows,{'expNPD21LS'})) ./ table2array(myMapelOptions.brTableRECIST(regRows,{'expN'}));
    weights = table2array(myMapelOptions.brTableRECIST(regRows,{'expN'}))/sum(table2array(myMapelOptions.brTableRECIST(regRows,{'expN'})));
    expSD = pred(:,1);
    expPD = pred(:,2);
    fitTable = table(expSD,expPD,resp);
    lmResults = fitlm(fitTable,'resp ~ expSD + expPD','Intercept',false,'Weights',weights);
    % Predicted expNPD21LS numbers
    % pred*lmResults.Coefficients{:,'Estimate'}.*table2array(myMapelOptions.brTableRECIST(regRows,{'expN'}))
    % Observed
    % table2array(myMapelOptions.brTableRECIST(regRows,{'expNPD21LS'}))
else
    lmResults = {};
    warning(['Unable to run ',mfilename,'.  Returning an empty object.'])
end
end

