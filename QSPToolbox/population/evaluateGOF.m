function myVPop = evaluateGOF(myVPop)
% Evaluate the goodness-of-fit for a given
% axis weight solution (VPop).  This has been implemented
% as a stand-alone function rather than a method.
%
% ARGUMENTS
% myVPop:         An instance of a VPop, VPopRECIST, or VPopRECISTnoBin object.  Populating the tables
%                 with experiment and simulated data should be performed
%                 prior to calling this function:
%                  mnSDTable
%                  binTable
%                  distTable
%                  distTable2D
%                  brTableRECIST (only if a VPopRECIST or VPopRECISTnoBin object)									   
%
% RETURNS
% myVpop:         A VPop is returned, with the fields populated:
%                  gofMn
%                  gofSD
%                  gofBin
%                  gofDist
%                  gofDist2D
%                  gof
%

continueFlag = false;
if nargin > 1
    warning(['Too many input arguments provided to,',mfilename,'.  Requires: myVPop.'])
    continueFlag = false;
elseif nargin > 0
    continueFlag = true;
else
    continueFlag = false;
    warning(['Insufficient input arguments provided to,',mfilename,'.  Requires: myVPop.'])
end

if continueFlag
    if ~ismember(class(myVPop),{'VPop','VPopRECIST','VPopRECISTnoBin'})
        continueFlag = false;
        warning(['Input VPop not recognized in call to ',mfilename,'.'])
    end       
end

if continueFlag
    % First, update the predicted values
    myVPop = myVPop.addPredTableVals();

    myMnSDTable = myVPop.mnSDTable;
    myBinTable = myVPop.binTable;
    myDistTable = myVPop.distTable;
	if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPopRECISTnoBin')
		myBRTableRECIST = myVPop.brTableRECIST;
		myRTableRECIST = myVPop.rTableRECIST;        
	end
	if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPopRECISTnoBin')
		myDistTable2D = myVPop.distTable2D;	
	else
		% Not implemented yet
		myDistTable2D = '';
	end
    

    if ~isempty(myMnSDTable)
        % Evaluate mn sd first
        meanResids = myMnSDTable{:,'expMean'} - myMnSDTable{:,'predMean'};
        expN = myMnSDTable{:,'expN'};
        predN = myMnSDTable{:,'predN'};
        expSD = myMnSDTable{:,'expSD'};
        predSD = myMnSDTable{:,'predSD'};

        % First we apply the t-test
        denominator = ( ((predN-1).*predSD.^2 + (expN-1).*expSD.^2)./(predN+expN-2) ).^.5;
        denominator = denominator.*((predN + expN)./(predN.*expN)).^.5;
        deltaMean = myMnSDTable{:,'predMean'} - myMnSDTable{:,'expMean'};
        tStat = deltaMean./denominator;
        % We want the two-sided test
        meanPvals = tcdf(-1*abs(tStat),predN+expN-2);
        meanPvals = 2 * meanPvals;

        % Next we apply the F-test
        sdPvals = fcdf((predSD./expSD).^2, predN-1, expN-1);

        % We want the two-sided test result
        sdPvals = 2 * min([sdPvals, 1-sdPvals],[],2);
        myVPop.gofMn = meanPvals;
        myVPop.gofSD = sdPvals;
    else
        myVPop.gofMn = [];
        myVPop.gofSD = [];
    end
    if ~isempty(myBinTable)
        expN = myBinTable{:,'expN'};
        predN = myBinTable{:,'predN'};
        expBinCounts = round(repmat(expN,1,4) .* [myBinTable{:, 'expBin1'}, myBinTable{:, 'expBin2'}, myBinTable{:, 'expBin3'}, myBinTable{:, 'expBin4'}]);
        predBinCounts = round(repmat(predN,1,4) .* [myBinTable{:, 'predBin1'}, myBinTable{:, 'predBin2'}, myBinTable{:, 'predBin3'}, myBinTable{:, 'predBin4'}]);
        [nStatRows, nBins] = size(expBinCounts);
        chiPvals = nan(nStatRows,1);
        for rowCounter = 1 : nStatRows 
            % The desired test is a 2x4 contigency table.  For sufficient smaple
            % sizes, chi^2 approximation is used.  For smaller sizes, an alternate
            % approximation to the exact test is employed.
            chiPval = contingency24(expBinCounts(rowCounter,:), predBinCounts(rowCounter,:),myVPop.exactFlag);
            % A 3rd party myfisher24() is used for the small sample size table.
            chiPvals(rowCounter) = chiPval;
        end
        % Write these pvals to the VPop before returning it.
        myVPop.gofBin = chiPvals;
    else
        myVPop.gofBin = [];
    end
    if ~isempty(myDistTable)
        [nStatRows, ~] = size(myDistTable);
        ksPvals = nan(nStatRows,1);
		predW = myDistTable.('predProbs');
		predInd = myDistTable.('simCombinedIndices');
		expInd = myDistTable.('expCombinedIndices');
        expN = myDistTable.('expN');
        predN = myDistTable.('predN');		
        for rowCounter = 1 : nStatRows 
            % Custom KS test for weighted samples
            expW = 1./(expN(rowCounter)*ones(1,expN(rowCounter)));
            expSample = myDistTable.('expSample'){rowCounter};
            predSample = myDistTable.('predSample'){rowCounter};
            SC = myDistTable.('combinedPoints'){rowCounter};            
            ksPval = weightedKSPreGrid(SC, expInd{rowCounter}, predInd{rowCounter}, expW, predW{rowCounter}, expN(rowCounter), predN(rowCounter));
            ksPvals(rowCounter) = ksPval;
        end		
        % Write these pvals to the VPop before returning it.
        myVPop.gofDist = ksPvals;
    else
        myVPop.gofDist = [];
	end
    if ~isempty(myDistTable2D)
        [nStatRows, ~] = size(myDistTable2D);
        ksPvals2D = nan(nStatRows,1);
		predW = myDistTable2D.('predProbs');
		predInd = myDistTable2D.('predIndices');
        expN = myDistTable2D.('expN');
        predN = myDistTable2D.('predN');
		combineQuadrants = myDistTable2D.('combinedQuadrants');
        for rowCounter = 1 : nStatRows 
			expW = 1./(expN(rowCounter)*ones(1,expN(rowCounter)));
            % Custom KS test for weighted samples
			quadrantCounts1c = combineQuadrants{rowCounter,1};
			quadrantCounts1s = combineQuadrants{rowCounter,2};
			quadrantCounts2c = combineQuadrants{rowCounter,3};
			quadrantCounts2s = combineQuadrants{rowCounter,4};
			ksPval = weightedKS2D(quadrantCounts1c, quadrantCounts1s, quadrantCounts2c, quadrantCounts2s, expW, predW{rowCounter}, expN(rowCounter), predN(rowCounter));
            ksPvals2D(rowCounter) = ksPval;
        end		
        % Write these pvals to the VPop before returning it.
        myVPop.gofDist2D = ksPvals2D;
    elseif isa(myVPop,'VPopRECIST') || isa(myVPop,'VPopRECISTnoBin')
        myVPop.gofDist2D = [];
	end	
	if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPopRECISTnoBin')
		if ~isempty(myBRTableRECIST)
			expN = myBRTableRECIST{:,'expN'};
			predN = myBRTableRECIST{:,'predN'};
			expBinCounts = round(repmat(expN,1,4) .* [myBRTableRECIST{:, 'expCR'}, myBRTableRECIST{:, 'expPR'}, myBRTableRECIST{:, 'expSD'}, myBRTableRECIST{:, 'expPD'}]);
			predBinCounts = round(repmat(predN,1,4) .* [myBRTableRECIST{:, 'predCR'}, myBRTableRECIST{:, 'predPR'}, myBRTableRECIST{:, 'predSD'}, myBRTableRECIST{:, 'predPD'}]);
			[nStatRows, nBins] = size(expBinCounts);
			chiPvals = nan(nStatRows,1);
			for rowCounter = 1 : nStatRows 
				% The desired test is a 2x4 contigency table.  For sufficient smaple
				% sizes, chi^2 approximation is used.  For smaller sizes, an alternate
				% approximation to the exact test is employed.
				chiPval = contingency24(expBinCounts(rowCounter,:), predBinCounts(rowCounter,:),myVPop.exactFlag);
				% A 3rd party myfisher24() is used for the small sample size table.
				chiPvals(rowCounter) = chiPval;
			end
			% Write these pvals to the VPop before returning it.
			myVPop.gofBR = chiPvals;
		else
			myVPop.gofBR = [];
        end
		if ~isempty(myRTableRECIST)
			expN = myRTableRECIST{:,'expN'};
			predN = myRTableRECIST{:,'predN'};
			expBinCounts = round(repmat(expN,1,4) .* [myRTableRECIST{:, 'expCR'}, myRTableRECIST{:, 'expPR'}, myRTableRECIST{:, 'expSD'}, myRTableRECIST{:, 'expPD'}]);
			predBinCounts = round(repmat(predN,1,4) .* [myRTableRECIST{:, 'predCR'}, myRTableRECIST{:, 'predPR'}, myRTableRECIST{:, 'predSD'}, myRTableRECIST{:, 'predPD'}]);
			[nStatRows, nBins] = size(expBinCounts);
			chiPvals = nan(nStatRows,1);
			for rowCounter = 1 : nStatRows 
				% The desired test is a 2x4 contigency table.  For sufficient smaple
				% sizes, chi^2 approximation is used.  For smaller sizes, an alternate
				% approximation to the exact test is employed.
				chiPval = contingency24(expBinCounts(rowCounter,:), predBinCounts(rowCounter,:),myVPop.exactFlag);
				% A 3rd party myfisher24() is used for the small sample size table.
				chiPvals(rowCounter) = chiPval;
			end
			% Write these pvals to the VPop before returning it.
			myVPop.gofR = chiPvals;
		else
			myVPop.gofR = [];
		end         
	end			   
    if isa(myVPop,'VPopRECIST') || isa(myVPop,'VPopRECISTnoBin')
		myVPop = compositeGOFRECIST(myVPop);
	else
		myVPop = compositeGOF(myVPop);
	end
else
    warning(['Unable to run ',mfilename,'.  Returning input.'])
end    

end