function interventionIDs = getInterventionIDs(myWorksheet)
% Get VP IDs from a model
% ARGUMENTS
% myWorksheet: a worksheet
%
% RETURNS
% interventionIDs: an array of cells with the intervention names

interventions = myWorksheet.interventions;
[nInterventions, dummy] = size(interventions);
interventionIDs = cell(1,nInterventions);
for theInterventionCounter = 1 : nInterventions
    curIntervention = interventions{theInterventionCounter};
    interventionIDs{theInterventionCounter} = curIntervention.('ID');
end
        
        
end