function simStruct = convertSimData(simData)
% Take a simResult and convert it to a format just with 2 fields:
% Name and Data.  Originally, I was going to use tables but ran into
% maximum name length issue.
%
% ARGUMENTS
% simData
%
% RETURNS
% simStruct
%
colnames = ['time';simData.DataNames];
datavalues = [simData.Time,simData.Data];
simStruct.Data = datavalues;
simStruct.Names = transpose(colnames);
end