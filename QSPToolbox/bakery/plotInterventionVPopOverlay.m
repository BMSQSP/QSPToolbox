function plotHandleVector = plotInterventionVPopOverlay(myWorksheet, myVPop, myPlotOptions, plotExisting, myPlotColor, myPlotLine)
% This function plots results for a single variable stored in a worksheet.
% The returned plots show the 5th, 50th, and 95th percentile.
%
% ARGUMENTS
% myWorksheet:      (required) Worksheet with simulation results
% myVPop:           (required) Virtual population with prevalence weight
%                   solution
% myPlotOptions:    Options structure to adjust the plotting
% plotExisting:     Whether to plot over the existing figure (true, false)
% myPlotColor:      Color scheme for plotting ('blue', 'lgrey', 'dgrey')
% myPlotLine:       MATLAB linestyle for plotting 
%
% RETURNS
% plotHandle
%

% Check number of inputs.
% We define fail behavior as warn rather than error here
% This just seems a more reasonable way to fail if we can generate plots
% rather than interrupting higher level execution with an error.
plotHandleVector = [];
flagExpData = false;

if nargin > 6
    flagContinue = false;
    warning(['Too many input arguments for ',mfilename,', require: myWorksheet, myVPop; optionally myPlotOptions, plotExisting, myPlotColor. Exiting.'])
elseif nargin > 5
    flagContinue = true;
elseif nargin > 4
    myPlotLine = '-';
    flagContinue = true;    
elseif nargin > 3
    myPlotLine = '-';    
    flagContinue = true;
    myPlotColor = 'blue';
elseif nargin > 2
    myPlotLine = '-';    
    flagContinue = true;
    plotExisting = false;
    myPlotColor = 'blue'    
else
    flagContinue = false;
    warning(['Insufficient input arguments for ',mfilename,' require: myWorksheet, myVPop; optionally myPlotOptions, plotExisting, myPlotColor. Exiting.'])  
end

if flagContinue
    myPlotColor = lower(myPlotColor);
    if sum(ismember({'blue','dgrey','lgrey'},myPlotColor)) ~= 1
        flagContinue = false;
        warning(['Invalid myPlotColor for ',mfilename,'.  Valid options are: blue, dgrey, lgrey.  Exiting.']);
    end                   
end

if flagContinue
    myPlotColor = lower(myPlotColor);
    if sum(islogical(plotExisting)) ~= 1
        flagContinue = false;
        warning(['Invalid plotExisting for ',mfilename,', exiting.']);
    end                   
end

if flagContinue
    if ~strcmp(class(myPlotOptions),'plotOptions')
        flagContinue = false;
        warning(['Invalid plotOptions object for ',mfilename,', exiting.'])
    end        
    if ~strcmp(class(myVPop), 'VPop')
        flagContinue = false;
        warning(['Invalid VPop for ',mfilename,'.'])
    end             
end

if flagContinue    
    if length(myVPop.pws) < 1
        flagContinue = false;
        warning(['VPop contains no prevalence weight results for ',mfilename,', exiting.'])
    end            
end


if flagContinue
    if length(myPlotOptions.vpIDs) < 1
        myPlotOptions.vpIDs = getVPIDs(myWorksheet);
    end
end

if flagContinue
    vpIDs = myPlotOptions.vpIDs;
    allVPIDs = getVPIDs(myWorksheet);
    if (sum(ismember(allVPIDs,vpIDs))<length(allVPIDs))
        warning(['All VPs are considered for the purpose of plotting the results of the VPop in ',mfilename,'. Cannot proceed with a subset, exiting.'])
        flagContinue = false;
    end
end

if flagContinue
    if ((length(myPlotOptions.expDataID) > 0) || (length(myPlotOptions.expDataTimeVar) > 0) || (length(myPlotOptions.expDataYVar) > 0))
        flagExpData = true;
    end
end

if flagContinue
    flagContinue = myPlotOptions.verify(myWorksheet);
    if ~flagContinue
        warning(['The plotOptions failed verification for ',mfilename,', exiting.'])
    end
end

   
if flagContinue
    nVPs = length(allVPIDs);
    interventionID = myPlotOptions.interventionID;    
    varName = myPlotOptions.varName;
end

if flagExpData && flagContinue
    % Filter out nan's, this will cause issues 
    expDataID = myPlotOptions.expDataID;
    expData = getExpData(expDataID, myWorksheet);
    expDataTimeVar = myPlotOptions.expDataTimeVar;
    expDataYVar = myPlotOptions.expDataYVar;
    expDataTime = expData.(expDataTimeVar);
    expDataYVar = expData.(expDataYVar);
    the_indices = find((~isnan(expDataTime)) & (~isnan(expDataYVar)));
    expDataTime = expDataTime(the_indices);
    expDataYVar = expDataYVar(the_indices);
end

if flagContinue
    % We reconstruct data here rather than getting from VPop in case
    % new variables will be plotted that weren't in the calibration
    interventionIDs = getInterventionIDs(myWorksheet);
    wshInterventionIndex = find(ismember(interventionIDs,interventionID));
    [nTimePoints, ~] = size(myWorksheet.results{wshInterventionIndex, 1}.Data);
    dataToPlot = nan(nTimePoints,nVPs);
    % Time should be constant
    for vpCounter = 1 : nVPs
        curResult = myWorksheet.results{wshInterventionIndex, vpCounter};
        curTimeIndex = find(ismember(curResult.Names,'time'));
        curVarIndex = find(ismember(curResult.Names,varName));
        curTime = curResult.Data(:,curTimeIndex);
        curData = curResult.Data(:,curVarIndex);
        % Interpolate to get the time exactly right
        % interpolateValue = interp1(curTime,curVar,expTime,'linear');
        dataToPlot(:, vpCounter) = curData;
    end
    
    % Plot 5/50/95 prediction interval
    intervalsToPlot = nan(nTimePoints,3);
    for timePointCounter = 1 : nTimePoints
        curData = dataToPlot(timePointCounter,:);
        [curData, indices] = sort(curData);
        curPWs = myVPop.pws(indices);
        curPWs = cumsum(curPWs);
        % PWs may be very small, and may be zero especially with rounding.
        % so we'll use both "first" and "last" occurrence of the PW comsum
        % and take the median of intervening values as the point for
        % to inform the percentile extrapolation
        [temp, indicesFirst] = unique(curPWs,'first');
        [curPWs, indicesLast] = unique(curPWs,'last');
        nUnique = length(indicesFirst);
        filteredData = nan(1, nUnique);
        for dataPointCounter = 1 : nUnique
            currentValues = curData(indicesFirst(dataPointCounter):indicesLast(dataPointCounter));
            filteredData(dataPointCounter) = median(currentValues);
        end
        intervalsToPlot(timePointCounter,:) = interp1(curPWs,filteredData,[0.05,0.50,0.95],'linear');
    end
    
    if ~plotExisting
        figure;
    end
    hold on
    if strcmp('blue',myPlotColor)
        color1 = 'b';
    elseif strcmp('dgrey',myPlotColor)
        color1 = [0.4,0.4,0.4];
    elseif strcmp('lgrey',myPlotColor)
        color1 = [0.7,0.7,0.7];
    end
    plotHandleVector(1) = fill([curTime-myPlotOptions.xShiftSim;flipud(curTime-myPlotOptions.xShiftSim)],[intervalsToPlot(:,1);flipud(intervalsToPlot(:,2))],color1,'edgecolor','none');
    set(plotHandleVector(1),'FaceAlpha',0.5);
    plotHandleVector(3) = fill([curTime-myPlotOptions.xShiftSim;flipud(curTime-myPlotOptions.xShiftSim)],[intervalsToPlot(:,2);flipud(intervalsToPlot(:,3))],color1,'edgecolor','none');
    set(plotHandleVector(3),'FaceAlpha',0.5); 
    plotHandleVector(1) = plot(curTime-myPlotOptions.xShiftSim, intervalsToPlot(:,1),'Color',color1,'LineWidth',4);
    plotHandleVector(2) = plot(curTime-myPlotOptions.xShiftSim, intervalsToPlot(:,2),myPlotLine,'Color',color1,'LineWidth',4);
    plotHandleVector(3) = plot(curTime-myPlotOptions.xShiftSim, intervalsToPlot(:,3),'Color',color1,'LineWidth',4);
    
    set(plotHandleVector(1),'linewidth',1)
    set(plotHandleVector(3),'linewidth',1)
    if (strcmp(myPlotOptions.scale,'xlogylin') || strcmp(myPlotOptions.scale,'xlogylog'))
        set(gca,'xscale','log'); % gca
    end
    if (strcmp(myPlotOptions.scale,'xlogylog') || strcmp(myPlotOptions.scale,'xlinylog'))
        set(gca,'yscale','log'); % gca
    end
    if flagExpData
        plotHandleVector(4) = plot(expDataTime, expDataYVar, 'k.','MarkerSize',20); 
        if (strcmp(myPlotOptions.scale,'xlogylin') || strcmp(myPlotOptions.scale,'xlogylog'))
            set(gca,'xscale','log'); % gca
        end
        if (strcmp(myPlotOptions.scale,'xlogylog') || strcmp(myPlotOptions.scale,'xlinylog'))
            set(gca,'yscale','log'); % gca
        end    
        %set(get(get(plotHandleVector(nVP+1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    if length(myPlotOptions.yLim) > 0
        ylim(myPlotOptions.yLim);
    end 
    if length(myPlotOptions.xLabelPretty) > 0
        xlabel(myPlotOptions.xLabelPretty);
    else
        xlabel(xVarName,'interpreter','none');
    end
    if length(myPlotOptions.yLabelPretty) > 0
        ylabel(myPlotOptions.yLabelPretty);
    else
        ylabel(varName,'interpreter','none');
    end    
    if length(myPlotOptions.xLim) > 0
        xlim(myPlotOptions.xLim);
    end    
    if length(myPlotOptions.yLim) > 0
        ylim(myPlotOptions.yLim);
    end     
    grid on
    set(gca,'fontsize', 18);
    if myPlotOptions.flagSave
        if length(myPlotOptions.fileName) > 0
            % saveas(gca,[myPlotOptions.fileName,'.tif']);
            print([myPlotOptions.fileName,'.tif'],'-dtiff','-r300');
        else
            theDate = date;
            formatOut = 'yymmdd';
            theDate = datestr(theDate,formatOut);            
            print([interventionID,'_',varName,'_vpop_',theDate,'.tif'],'-dtiff','-r300');
        end
    end
end
end