function plotHandle = plotBinVPop(myVPop, myPlotOptions)
% This function takes the bin data in a VPop distribution table
% to plot the experimental data side-by-side with the population
% to facilitate diagnosing fitting during
% the population calibration/prevalence weighting process.
%
% ARGUMENTS
% myVPop:           (required) A VPopRECIST or VPopRECISTnoBin object.
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
    if ~(ismember(class(myVPop), {'VPopRECIST','VPopRECISTnoBin'}))
        flagContinue = false;
        warning(['Invalid VPopRECIST for ',mfilename,'.'])
    end  
    if ~strcmp(class(myPlotOptions),'plotOptions')
        flagContinue = false;
        warning(['Invalid plotOptions for ',mfilename,'.'])
    end        
end

if flagContinue
    if ~(strcmp(class(myVPop.brTableRECIST), 'table'))
        flagContinue = false;
        warning(['Invalid VPop.brTableRECIST for ',mfilename,'.'])
    end  
end

if (flagContinue)
    [myNPlots, ~] = size(myVPop.brTableRECIST);
    if myNPlots < 1
        flagContinue = false;
        warning(['Invalid VPop.brTableRECIST for ',mfilename,'.'])        
    end
end


if (flagContinue)
    myTable = myVPop.brTableRECIST;
    
    for rowCounter = 1 : myNPlots
        plotNames{rowCounter} = {myTable{rowCounter,'interventionID'}{1},myTable{rowCounter,'elementID'}{1},num2str(myTable{rowCounter,'time'})};
    end
    
end

if (flagContinue)
    nPlotsVer = floor(sqrt(myNPlots));
    nPlotsHor = ceil(myNPlots / nPlotsVer);
    plotHandle = figure;
    for plotCounter = 1 : myNPlots
        subplot(nPlotsVer, nPlotsHor, plotCounter);
        expBins = [myTable{plotCounter,'expCR'},myTable{plotCounter,'expPR'},myTable{plotCounter,'expSD'},myTable{plotCounter,'expPD'}];
        predBins = [myTable{plotCounter,'predCR'},myTable{plotCounter,'predPR'},myTable{plotCounter,'predSD'},myTable{plotCounter,'predPD'}];
        binLabels = {'CR','PR','SD','PD'};
        barplot = bar([expBins;predBins]');
        set(gca,'XTickLabel',binLabels);
        xtickangle(45);
        set(barplot(1),'facecolor','w') 
        set(barplot(2),'facecolor',[.2,.2,.2]) 
        legend({'Trials','VPop'})
        title(gca,plotNames{plotCounter},'interpreter','none','FontWeight','Normal');
        set(gca,'box','on');
        set(gca,'fontsize', 10);
    end

    if myPlotOptions.flagSave
        theDate = date;
        formatOut = 'yymmdd';
        theDate = datestr(theDate,formatOut);
        print(['VPopBR_',theDate,'.tif'],'-dtiff','-r300');
    end
end

end