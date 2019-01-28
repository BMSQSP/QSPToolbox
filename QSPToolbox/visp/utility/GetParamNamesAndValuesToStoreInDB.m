function [paramValues,paramNames] = GetParamNamesAndValuesToStoreInDB(modelValues,modelNames,indices)

paramValues = modelValues(:,:,indices.varyingParam); 
paramNames  = modelNames(indices.varyingParam);

end


