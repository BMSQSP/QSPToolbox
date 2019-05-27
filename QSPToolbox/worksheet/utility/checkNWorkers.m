function mySimulateOptions = checkNWorkers(mySimulateOptions)
% This is a "utility function" that should be called directly.
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
if isnan(mySimulateOptions.nWorkers)
	mySettings = parallel.Settings;
	mySchedulerComponent = [mySimulateOptions.clusterID,'SchedulerComponent'];
	for checkCounter = 1 : length(mySettings.SchedulerComponents)
		if isequal(mySettings.SchedulerComponents(checkCounter).Name,mySchedulerComponent)
            % 'Use Default' is sometimes given
            if ~isnumeric(mySettings.SchedulerComponents(checkCounter).NumWorkers)
                mySimulateOptions.nWorkers = mySettings.Pool.PreferredNumWorkers;
            else
                mySimulateOptions.nWorkers = mySettings.SchedulerComponents(checkCounter).NumWorkers;
            end
        end
	end	
end
end
