function plotHandle = plotBinVPop(myVPop, myPlotOptions)
% This function takes the bin data in a VPop distribution table
% to plot the experimental data side-by-side with the population
% to facilitate diagnosing fitting during
% the population calibration/prevalence weighting process.
%
% ARGUMENTS
% myVPop:           (required) A VPop, VPopRECIST, or VPopRECISTnoBin.  Note  
%                   that both the experimental and simulation result fields should  
%                   be populated
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
    if ~(strcmp(class(myVPop.binTable), 'table'))
        flagContinue = false;
        warning(['Invalid VPop.binTable for ',mfilename,'.'])
    end  
end

if (flagContinue)
    [myNPlots, ~] = size(myVPop.binTable);
    if myNPlots < 1
        flagContinue = false;
        warning(['Invalid VPop.binTable for ',mfilename,'.'])        
    end
end


if (flagContinue)
    myTable = myVPop.binTable;
    
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
        expBins = [myTable{plotCounter,'expBin1'},myTable{plotCounter,'expBin2'},myTable{plotCounter,'expBin3'},myTable{plotCounter,'expBin4'}];
        predBins = [myTable{plotCounter,'predBin1'},myTable{plotCounter,'predBin2'},myTable{plotCounter,'predBin3'},myTable{plotCounter,'predBin4'}];
        binEdges = [myTable{plotCounter,'binEdge1'},myTable{plotCounter,'binEdge2'},myTable{plotCounter,'binEdge3'}];
        binLabels = {['< ',num2str(binEdges(1))],[num2str(binEdges(1)),' - ',num2str(binEdges(2))],[num2str(binEdges(2)),' - ',num2str(binEdges(3))],['>= ',num2str(binEdges(3))]};
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
        print(['VPopBin_',theDate,'.tif'],'-dtiff','-r300');
    end
end

end