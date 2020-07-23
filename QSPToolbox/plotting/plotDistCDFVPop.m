function plotHandle = plotDistCDFVPop(myVPop, myPlotOptions)
% This function takes the distribution data in a VPop distribution table
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
    if ~(strcmp(class(myVPop.distTable), 'table'))
        flagContinue = false;
        warning(['Invalid VPop.distTable for ',mfilename,'.'])
    end  
end

if (flagContinue)
    [myNPlots, ~] = size(myVPop.distTable);
    if myNPlots < 1
        flagContinue = false;
        warning(['Invalid VPop.distTable for ',mfilename,'.'])        
    end
end


if (flagContinue)
    myTable = myVPop.distTable;
    
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
        expN = myTable{plotCounter,'expN'};
        predN = myTable{plotCounter,'predN'};
        expW = 1./(expN*ones(1,expN));
        predW = myTable{plotCounter,'predProbs'}{1};
        expSample = myTable{plotCounter,'expSample'}{1};
        predSample = myTable{plotCounter,'predSample'}{1};
		SC = myTable{plotCounter,'combinedPoints'}{1};
        expInd = myTable{plotCounter,'expCombinedIndices'}{1};;
        predInd = myTable{plotCounter,'simCombinedIndices'}{1};;		
		[CDF1, CDF2] = alignCDFsPreGrid(SC, expInd, predInd, expW, predW);
        plot(SC,CDF1,'k-',SC,CDF2,'r-','LineWidth',3);
        ylim([0 1]);
        xlim([min(SC) max(SC)]);
        grid on
        title(gca,plotNames{plotCounter},'interpreter','none','FontWeight','Normal');
        set(gca,'box','on');
        set(gca,'fontsize', 10);
        
    end

    if myPlotOptions.flagSave
        theDate = date;
        formatOut = 'yymmdd';
        theDate = datestr(theDate,formatOut);
        print(['VPopDistCDF_',theDate,'.tif'],'-dtiff','-r300');
    end
end

end