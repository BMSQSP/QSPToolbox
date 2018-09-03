function interventionSubElements = extractInterventionTypeElements(intervention, type)
% Get information from an intervention
% ARGUMENTS
% intervention: an intervention structure
% elements of type to return: e.g. 'VARIANT', 'DOSE'
%
% RETURNS
% interventionSubElements: a cell array with the selection
%
interventionSubElements = cell(0,1);
[nModifiers, nCols] = size(intervention.('definition'));
for modifierCounter = 1 : nModifiers
    specType = intervention.('definition'){modifierCounter,1};
    specName = intervention.('definition'){modifierCounter,2};
    if strcmp(specType, type)
        interventionSubElements{length(interventionSubElements) + 1} = specName;
    end
end
end