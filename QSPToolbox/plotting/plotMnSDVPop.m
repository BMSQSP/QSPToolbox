function plotHandle = plotMnSDVPop(myVPop, myPlotOptions)
% This function takes the mean & SD matching data in a VPop MnSD table
% to plot the experimental data side-by-side with the population
% to facilitate diagnosing fitting during
% the population calibration/prevalence weighting process.
%
% ARGUMENTS
% myVPop:           (required) A VPop.  Note that both the 
%                   experimental and simulation result fields should be 
%                   populated
% myPlotOptions:    (optional) A plotOptions structure.
%                   Note that not all arguments are used.
%
% RETURNS
% plotHandle
%
plotHandle = [];
flagContinue = false;

if nargin > 2
    flagContinue = false;
    warning(['Too many input arguments to ',mfilename,'. Required: myVPop.'])
elseif nargin > 1
    flagContinue = true;    
elseif nargin > 0
    myPlotOptions = plotOptions;    
    flagContinue = true;
else
    flagContinue = false;
    warning(['Insufficient input arguments to ',mfilename,'. Required: myVPop.'])
end

if flagContinue
    if ~(ismember(class(myVPop), {'VPop','VPopRECIST','VPopRECISTnoBin'}))
        flagContinue = false;
        warning(['Invalid VPop for ',mfilename,'.'])
    end  
    if ~strcmp(class(myPlotOptions),'plotOptions')
        flagContinue = false;
        warning(['Invalid plotOptions for ',mfilename,'.'])
    end        
end

if flagContinue
    if ~(strcmp(class(myVPop.mnSDTable), 'table'))
        flagContinue = false;
        warning(['Invalid VPop.mnSDTable for ',mfilename,'.'])
    end  
end

if (flagContinue)
    [myNPlots, ~] = size(myVPop.mnSDTable);
    if myNPlots < 1
        flagContinue = false;
        warning(['Invalid VPop.mnSDTable for ',mfilename,'.'])        
    end
end

if (flagContinue)
    myMnValExp = myVPop.mnSDTable{:,'expMean'};
    myUBValExp = myVPop.mnSDTable{:,'expSD'};
    myLBValExp = myVPop.mnSDTable{:,'expSD'};
    
    myMnValVPop = myVPop.mnSDTable{:,'predMean'};
    myUBValVPop = myVPop.mnSDTable{:,'predSD'};
    myLBValVPop = myVPop.mnSDTable{:,'predSD'};
    
    for rowCounter = 1 : myNPlots
        plotNames{rowCounter} = {myVPop.mnSDTable{rowCounter,'interventionID'}{1},myVPop.mnSDTable{rowCounter,'elementID'}{1},num2str(myVPop.mnSDTable{rowCounter,'time'})};
    end
    
end

if (flagContinue)
    nPlotsVer = floor(sqrt(myNPlots));
    nPlotsHor = ceil(myNPlots / nPlotsVer);
    plotHandle = figure;
    for plotCounter = 1 : myNPlots
        subplot(nPlotsVer, nPlotsHor, plotCounter);
        fill([-2,2,2,-2],[myMnValExp(plotCounter)+myUBValExp(plotCounter),myMnValExp(plotCounter)+myUBValExp(plotCounter),myMnValExp(plotCounter)-myLBValExp(plotCounter),myMnValExp(plotCounter)-myLBValExp(plotCounter)], [0.75 0.75 0.75])
        hold on
        errorbar(-1,myMnValExp(plotCounter),myLBValExp(plotCounter),myUBValExp(plotCounter),'ko-','MarkerFaceColor','k','LineWidth',3);
        hold on
        errorbar(1,myMnValVPop(plotCounter),myLBValVPop(plotCounter),myUBValVPop(plotCounter),'ro-','MarkerFaceColor','r','LineWidth',3);
        hold on
        xlim([-1.25 1.25]);
        
        set(gca,'TickLabelInterpreter', 'none');
        set(gca, 'XTick', [-1 1]);
        set(gca, 'XTickLabel', {'Exp','VPop'});
        set(gca,'TickLabelInterpreter', 'none');
        grid on
        title(gca,plotNames{plotCounter},'interpreter','none','FontWeight','Normal');
        
        set(gca,'box','on');
        set(gca,'fontsize', 10);
        
    end

    if myPlotOptions.flagSave
        theDate = date;
        formatOut = 'yymmdd';
        theDate = datestr(theDate,formatOut);
        print(['VPopMnSD_',theDate,'.tif'],'-dtiff','-r300');
    end
end

end