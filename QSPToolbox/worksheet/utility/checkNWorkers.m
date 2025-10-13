function mySimulateOptions = checkNWorkers(mySimulateOptions)
% This is a "utility function" that should not be called directly.
% This is used to check the number of workers specified
% and if not, check how many are available.
%
% ARGUMENTS
%  mySimulateOptions:          a simulateOptions object
%
% RETURNS
%  mySimulateOptions:          a simulateOptions object with the
%                              number of workers checked
%


mySettings = parallel.Settings;
clusterCheck = parcluster(mySimulateOptions.clusterID);
if isnan(mySimulateOptions.nWorkers)
    % The purpose of this block is to assign to nWorkers the value of
    % maxNumWorkers, preferredNumWorkers, numWorkers... whatever field that is
    % available in the clusterCheck object
    if ismember(version("-release"), {'2024a'})  % changed from 2024a
        mySimulateOptions.nWorkers = clusterCheck.PreferredPoolNumWorkers;
    else
        if isnumeric(mySettings.Pool.PreferredNumWorkers)
            % First try getting the value from parallel.settings, but cap in
            % case a silly value was entered.
            if isa(clusterCheck, 'parallel.cluster.Local')
                mySimulateOptions.nWorkers = min(clusterCheck.NumWorkers, ...
                    mySettings.Pool.PreferredNumWorkers);
            elseif (strcmp(clusterCheck.Profile,'mps-slurm') | strcmp(clusterCheck.Profile,'mps-hybrid') )
                mySimulateOptions.nWorkers = min(clusterCheck.NumWorkers, ...
                    mySettings.Pool.PreferredNumWorkers);
            else
                mySimulateOptions.nWorkers = min(clusterCheck.MaxNumWorkers, ...
                    mySettings.Pool.PreferredNumWorkers);
            end
        else
            mySimulateOptions.nWorkers = clusterCheck.NumWorkers;
        end
    end
else
    if ismember(version("-release"), {'2024a'})
        mySimulateOptions.nWorkers = min(mySimulateOptions.nWorkers, clusterCheck.PreferredPoolNumWorkers);
    else
        %     Otherwise, check that the number of workers specified
        %     doesn't exceed any pre-set maximums.
        if isnumeric(mySettings.Pool.PreferredNumWorkers)
            % First try getting the value from parallel.settings, but cap in
            % case a silly value was entered.
            if isa(clusterCheck, 'parallel.cluster.Local')
                mySimulateOptions.nWorkers = min(mySimulateOptions.nWorkers, min(clusterCheck.NumWorkers, mySettings.Pool.PreferredNumWorkers));
            else
                mySimulateOptions.nWorkers = min(mySimulateOptions.nWorkers, min(clusterCheck.MaxNumWorkers, mySettings.Pool.PreferredNumWorkers));
            end
        else
            if isa(clusterCheck, 'parallel.cluster.Local')
                mySimulateOptions.nWorkers = min(mySimulateOptions.nWorkers, clusterCheck.NumWorkers);
            else
                mySimulateOptions.nWorkers = min(mySimulateOptions.nWorkers, clusterCheck.MaxNumWorkers);
            end
        end
    end
end
