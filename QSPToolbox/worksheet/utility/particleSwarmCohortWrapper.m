function [x,fval,exitFlag,swarm] = particleSwarmCohortWrapper(fun, nvars, lb, ub, userOptions)
% This is a wrapper to call particle swarm and keep the last returned
% swarm, which isn't available in the userOptions.

swarm = [];
userOptions.OutputFcn = @outputFunction;
[x,fval,exitFlag,output] = particleswarm(fun, nvars, lb, ub, userOptions);

    function stop = outputFunction(optimValues, state)
        stop = false;
        swarm = optimValues.swarm;
    end
end