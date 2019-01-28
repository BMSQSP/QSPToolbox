function simBioModel = UpdateSimBioModel(modelElements,simBioModel,modelValues,indices)
  for counter = 1:length(indices.modelElemToSet)
    paramName      = modelElements{indices.modelElemToSet(counter),1};
    paramType      = modelElements{indices.modelElemToSet(counter),2};
    paramValue     = modelValues(1,1,indices.constantParam(counter));
    myObject       = sbioselect(simBioModel,'Type',paramType,'Name',paramName);
    set(myObject,'Value',paramValue); 
  end
end

