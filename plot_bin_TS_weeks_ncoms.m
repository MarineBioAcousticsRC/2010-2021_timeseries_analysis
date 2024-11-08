clearvars

minClicksSp =  {'Boat';'Gg';'Kspp';'Md';'Me';'Noise';'Pm';'Sonar';'UD';'UD_3p';...
    'UD_LF';'Zc';'UD_all';'SS'};

minClicks = 10; % minimum number of clicks required to count a bin as positive.
myThreshVals = 0.5; % Minimum label confidence
fullReRun1 = 1;
siteVec = {'GC','MC','DT','DC','MP'};
overallStart = datenum(2010,1,1);
overallEnd = datenum(2021,09,20);
plotsOn = 1;

for iSite = 1:5
    nameCode = siteVec{iSite};

    %%
    if strcmp(nameCode,'GC')
        outDir = 'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Timeseries\GC\';
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        outStr = 'GOM_GC_bin';
        outName = fullfile(outDir,'GOM_GC_bin_TS_noNorm_delphTypes.mat');
        siteName = 'Green Canyon';
        siteCode = 'GC';
        myEffortFile = ('P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Effort\GC_Effort.mat');
        binDirList = {'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Labels\GC_ID'};
        concatzIDtimes = fullfile(outDir,'GC_concatZID.mat');

    elseif strcmp(nameCode,'MC')
        outDir = 'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Timeseries\MC\';
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        outStr = 'GOM_MC_bin';
        outName = fullfile(outDir,'GOM_MC_bin_TS_noNorm_delphTypes.mat');
        siteName = 'Mississippi Canyon';
        siteCode = 'MC';
        myEffortFile = ('P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Effort\MC_Effort.mat');
        binDirList = {'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Labels\MC_ID'};
        concatzIDtimes = fullfile(outDir,'MC_concatZID.mat');

    elseif strcmp(nameCode,'DT')
        outDir = 'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Timeseries\DT\';
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        outStr = 'GOM_DT_bin';
        outName = fullfile(outDir,'GOM_DT_bin_TS_noNorm_delphTypes.mat');
        siteName = 'Dry Tortugas';
        siteCode = 'GC';
        myEffortFile = ('P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Effort\DT_Effort.mat');
        binDirList = {'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Labels\DT_ID'};
        concatzIDtimes = fullfile(outDir,'DT_concatZID.mat');

    elseif strcmp(nameCode,'DC')
        outDir = 'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Timeseries\DC\';
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        outStr = 'GOM_DC_bin';
        outName = fullfile(outDir,'GOM_DC_bin_TS_noNorm_delphTypes.mat');
        siteName = 'DeSoto Canyon';
        siteCode = 'DC';
        myEffortFile = ('P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Effort\DC_Effort.mat');
        binDirList = {'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Labels\DC_ID'};
        concatzIDtimes = fullfile(outDir,'DC_concatZID.mat');

    elseif strcmp(nameCode,'MP')
        outDir = 'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Timeseries\MP\';
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        outStr = 'GOM_MP_bin';
        outName = fullfile(outDir,'GOM_MP_bin_TS_noNorm_delphTypes.mat');
        siteName = 'Main Pass';
        siteCode = 'MP';
        myEffortFile = ('P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Effort\MP_Effort.mat');
        binDirList = {'P:\doi_10_5061_dryad_9zw3r22n4__v20241022\GOM_2010-2021\Labels\MP_ID'};
        concatzIDtimes = fullfile(outDir,'MP_concatZID.mat');

    end
    %%
    cPath = {};
    cName = {};

    %% merge labels into one vector - this takes awhile.
    myDayBinsAll = [];
    myTSbinAll = [];
    myTSbinNoScalingAll = [];
    myTimeBinsAll = [];
    labelScore = [];
    %%
    if fullReRun1
        for iDep = 1:length(binDirList)
            zIDAllBin = [];
            thisBinDir = binDirList{iDep};

            dirList = dir(fullfile(thisBinDir,'\G*ID1.mat'));
            % labelDir = dir(fullfile(fileparts(thisBinDir),'*labels.mat'));
            labelDir =[];
            if isempty(labelDir)
                labelDir = dir(fullfile(strrep(thisBinDir,'_ID',''),'\G*labels.mat'));
            end

            if isempty(dirList)
                continue
            end


            fprintf('Beginning folder %s\n', thisBinDir)
            for iDir = 1:length(dirList)

                load(fullfile(thisBinDir,dirList(iDir).name))
                cPath = [cPath;classificationInfo.networkPath];
                cName = [cName;dirList(iDir).name];
                zIDAllBin = [zIDAllBin;zID];
                fprintf('     done with file %s\n', dirList(iDir).name)
            end

            % make vector of 5 min bins spanning the whole dataset
            myTimeBins = floor(min(zIDAllBin(:,1))):1/(24*12):ceil(max(zIDAllBin(:,1)));
            myTimeBinsAll = [myTimeBinsAll,myTimeBins];

            udIdx = (find(strcmp('UD',mySpID)));

            tpIdx = (find(strcmp('UD_3p',mySpID)));

            if ~isempty(tpIdx)
                % change all low probability ud3p to generic UD. 
                change3p = find(zIDAllBin(:,2)==tpIdx);
                zIDAllBin(change3p,2) = udIdx;
            end

            %countAdjust = bNClicks./bNClass;
            nSp = length(mySpID);
            if nSp >12
                1;
            end

            myTSbin = zeros(size(myTimeBins,2),nSp,length(myThreshVals));
            myTSbinNoScaling = zeros(size(myTimeBins,2),nSp,length(myThreshVals));
            fprintf('Computing time series\n')
            %[~,iTS,iCountAdj] = intersect(myTimeBins,bTimes(:,1));
            % Apply classification score threshold,
            for iSP = 1:nSp
                for iThresh = 1%:length(myThreshVals)
                    [I,~,tempTS] = histcounts(zIDAllBin((zIDAllBin(:,2)==iSP)&...
                        (zIDAllBin(:,3)>=myThreshVals(iThresh)),1),...
                        myTimeBins);
                    myTSbin(1:end-1,iSP,iThresh) = I;%.*countAdjust(iCountAdj));
                    myTSbinNoScaling(1:end-1,iSP,iThresh) = I;
                end
            end
            myTSbinAll = [myTSbinAll;myTSbin];
            myTSbinNoScalingAll = [myTSbinNoScalingAll;myTSbinNoScaling];

            myDayBins = floor(min(zIDAllBin(:,1))):1:ceil(max(zIDAllBin(:,1)));
            myDayBinsAll = [myDayBinsAll,myDayBins];
        end

        % save timeseries of bins from all files.
        save(concatzIDtimes,'myDayBinsAll','myTSbinNoScalingAll','myTSbinAll','labelScore',...
            'myTimeBinsAll','mySpID','nSp','-v7.3')
    else
        load(concatzIDtimes)

        if ~isfolder(outDir)
            mkdir(outDir)
        end
    end

    load(myEffortFile)
    dIdx = (find(strcmp('Dolphins',mySpID)));
    if ~isempty(dIdx)
        mySpID{dIdx} = 'UD';
    end
    [~,tsSort] = sort(myTimeBinsAll);
    myTimeBinsAllSort = myTimeBinsAll(tsSort);
    myTSbinAllSort = myTSbinAll(tsSort,:,:);
    myTSbinNoScalingAllSort = myTSbinNoScalingAll(tsSort,:,:);

    [~,tsDaySort] = sort(myDayBinsAll);
    myDayBinsAllSort = myDayBinsAll(tsDaySort);

    %% apply effort
    [effort,latLongs,errRates,depl] = gofmx_dates4(siteCode);
    earlyEffort = effort(:,1)<overallStart;
    effort(earlyEffort,:)=[];

    allEffort(allEffort(:,2)>1,2) = 1;
    allEffort(allEffort(:,2)<=0.5,2)= NaN;
    allEffort(:,1) = allEffort(:,1)+datenum([2000,0,0]);
    for iL=1:size(effort,1)-1

        outOfEffort = find(allEffort(:,1)>effort(iL,2)& allEffort(:,1)<effort(iL+1,1));
        if ~isempty(outOfEffort)
            allEffort(outOfEffort,2)=NaN;
        end
    end

    outOfEffortBefore = find(allEffort(:,1)<effort(1,1));
    allEffort(outOfEffortBefore,2)=NaN;
    outOfEffortAfter = find(allEffort(:,1)>effort(end,2));
    allEffort(outOfEffortAfter,2)=NaN;
    allEffort(~isnan(allEffort(:,2)),2) = 1;

    allEffort4Diel = [];
    s1 = 1;
    for iEff = 1:size(allEffort,1)-1
        myDiff = round(allEffort(iEff:iEff+1,2));
        if (myDiff(1)==1) && isnan(myDiff(2))
            allEffort4Diel(s1,2) = allEffort(iEff,1);
            s1 = s1+1;
        elseif isnan(myDiff(1)) && (myDiff(2)==1)
            allEffort4Diel(s1,1) = allEffort(iEff+1,1);
        end
    end
    save(outName,'allEffort','myDayBinsAllSort','myTSbinAllSort','myTimeBinsAllSort','mySpID',...
        'myTSbinNoScalingAllSort','myThreshVals')


    %% for slides
    nVals = 1;
    annualChangeMat = [];
    cvDensity = [];
    densityTimesAll = {};
    densityEstAll = {};
    densityCVAll = {};
    deseasTrendPlusMeanAll = {};
    allTScountsEffortAdjust = [];
    nktTkt = [];
    nkt = [];
    fpRateAll = [];
    fnRateAll = [];
    pDetAll = [];
    pVocalAll = [];
    grpSizeAll = [];
    startValAll = [];
    endValAll = [];
    bootSlopeAll = [];
    percentChange_10yrAll = [];

    for mySP = 1:nSp
        if plotsOn
            figure(300);clf
        end
        myWeekBinsAllSort = (overallStart:7:overallEnd);
        %         indata = overallStart:overallEnd;
        %         dvec   = datevec(indata);
        %         duniq  = unique(dvec(:, 1:2), 'rows');
        %         myWeekBinsAllSort = datenum(duniq(:,1), duniq(:,2), 1);
        weekMean = zeros(size(myWeekBinsAllSort,2),nVals);
        weekStd = zeros(size(myWeekBinsAllSort,2),nVals);
        for iB = 1:nVals

            minClickIdx = find(strcmp(mySpID{mySP},minClicksSp));
            if ~isempty(minClickIdx)
                minClickVal = minClicks;
            else
                minClickVal = 5;
            end
            [thisTScounts,~] = histc(myTimeBinsAllSort(myTSbinAllSort(:,mySP,1)>=minClickVal),myDayBinsAllSort);
            [~,goodCols,goodColsB] = intersect(mySpID,{'Pm','Zc','Me','Md','Kspp','Gg','UD','UD_LF'});
            [~,I] = unique(myDayBinsAllSort,'last');
            [intersectDates,iA,iC] = intersect(myDayBinsAllSort(I),allEffort(:,1));
            thisTScountsEffortAdjust = thisTScounts(I(iA))./allEffort(iC,2)';
            allTScountsEffortAdjust(:,mySP) = thisTScountsEffortAdjust;

            if mySP == length(mySpID)
                mergedPosDays = sum(sum(allTScountsEffortAdjust(:,goodCols),2)>0);
                mergedDays = sum(sum(allTScountsEffortAdjust(:,goodCols),2)>=0);
                mergedSpPosDays = [];
                for iGood = 1:length(goodCols)
                    mergedSpPosDays(1,goodColsB(iGood)) = sum(allTScountsEffortAdjust(:,goodCols(iGood))>0);
                end
            end

            [~,theseDays] = histc(intersectDates,myWeekBinsAllSort);
            weekEffort = [];
            for iBars = 1:length(myWeekBinsAllSort)

                dayInd = find(theseDays==iBars);
                if sum(~isnan(thisTScountsEffortAdjust(dayInd)))<5 ||...
                        length(dayInd)<=5
                    weekMean(iBars,iB) = NaN;
                    weekStd(iBars,iB) = NaN;
                    weekEffort(iBars,1) = NaN;
                    continue
                else
                    1;
                end
                weekMean(iBars,iB) = nanmean(thisTScountsEffortAdjust(dayInd));
                weekStd(iBars,iB) = nanstd(thisTScountsEffortAdjust(dayInd));
                weekEffort(iBars,1) = 1*12*24;%nanmean(allEffort(dayInd,2))*12*24; % note this is average bins per day, taken over the week.
                % presence should be computed the same way.
            end
            %             weekMean(end+1) = NaN; % fixing issue where last bin is zero not NaN - why?
            %             weekStd(end+1) = NaN;
            %             weekEffort(end+1,1) = NaN;
            if plotsOn
                subplot(nVals,1,iB)

                errorbar(myWeekBinsAllSort,weekMean(:,iB)*5,zeros(size(weekStd)),weekStd*5,'.k','CapSize',6);
                hold on

                % title(sprintf('%s: %s (bin labels)',siteName, mySpID{mySP}))
                t = annotation('textbox');
                t.String = siteCode;
                t.Position = [.86,.8,.05,.1];t.LineStyle = 'none';
                t.FontSize = 12; t.FontWeight = 'bold';

                xlabel('Date')
                ylabel('Weekly mean min/day')
                xlim([overallStart,overallEnd])

                tickList = get(gca,'xtick');
                thisYear = year(tickList)';
                thisMonth = month(tickList)';
                thisDay = ones(size(thisMonth));
                thisStartMonth = datenum([thisYear,thisMonth,thisDay]);
                set(gca,'xtick',thisStartMonth)
                datetick('x','mm/yy','keeplimits');grid on

                rH = get(gca,'ylim');
                bar(myWeekBinsAllSort(isnan(weekMean)),rH(2).*ones(size(weekMean(isnan(weekMean)))),...
                    1,'EdgeColor','none','FaceColor',[.76,.74,.99],'linestyle','none')
                if overallEnd-overallStart>500
                    myEdgeColor = 'b';
                else
                    myEdgeColor = 'k';
                end
                bar(myWeekBinsAllSort,weekMean(:,iB)*5,1,'Facecolor','b','EdgeColor',myEdgeColor)

                %end
                set(gca,'ylim',rH,'fontweight','bold')
                set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[2 2 10 5])

                set(gca,'layer','top')

                print(gcf,'-dpng','-r300',fullfile(outDir, [mySpID{mySP},'_',siteCode,'_noNorm_weekly_bin_Presence.png']))
                saveas(gcf,fullfile(outDir, [mySpID{mySP},'_',siteCode,'_noNorm_weekly_bin_Presence.fig']))
            end
        end

        uDates = overallStart:1:overallEnd;
        [sRise,sSet,~]=sunrise(latLongs(1,1),latLongs(1,2),0,0,uDates);

        %%
        if sum(strcmp(mySpID{mySP},{'Zc','Pm','Me','Md','Gg','UD','UD_LF','Kspp'}))>0

            [dayCount,I] = histc(myTimeBinsAllSort,[sRise;sSet]);
            [dayPosCount,I] = histc(myTimeBinsAllSort(myTSbinAllSort(:,mySP,1)>=minClickVal),[sRise;sSet]);

            dayBinsAll = (dayCount(1:2:end));
            nightBinsAll = (dayCount(2:2:end));
            daylightPerc = dayBinsAll./(dayBinsAll+nightBinsAll);
            badDaylight1 = find(daylightPerc<.4);
            badDaylight2 = find(daylightPerc>.6);
            daylightPerc(badDaylight1) = min(daylightPerc(badDaylight1-1:min(badDaylight1+1,length(daylightPerc))));
            daylightPerc(badDaylight2) = max(daylightPerc(badDaylight2-1:min(badDaylight2+1,length(daylightPerc))));
            [dayC, dayI] = histc(sRise,myWeekBinsAllSort);
            dayLightPercBin = [];
            for iBin = 1:length(myWeekBinsAllSort)
                dayLightPercBin(1,iBin) = nanmedian(daylightPerc(dayI==iBin));
            end
            dayPosSet = (dayPosCount(1:2:end));
            nightPosSet = (dayPosCount(2:2:end));

            dayPosSetNoNaN = dayPosSet(~isnan(dayBinsAll));
            nightPosSetNoNaN = nightPosSet(~isnan(nightBinsAll));

            nBins = length(dayPosSetNoNaN);
            dayRatio= [];
            nightRatio = [];
            for iBoot =1:100
                p = randperm(nBins,round(nBins*.25));
                dayRatio(iBoot) = sum(dayPosSetNoNaN(p))./(sum(dayPosSetNoNaN(p))+sum(nightPosSetNoNaN(p)));
                nightRatio(iBoot) = sum(nightPosSetNoNaN(p))./(sum(dayPosSetNoNaN(p))+sum(nightPosSetNoNaN(p)));
            end
            meanDayRatio(mySP,1) = nanmean(dayRatio(~isinf(dayRatio)));
            meanNightRatio(mySP,1) = nanmean(nightRatio(~isinf(dayRatio)));

            cvDayRatio(mySP,1) = nanstd(dayRatio(~isinf(dayRatio)))./meanNightRatio(mySP,1);
            cvNightRatio(mySP,1) = nanstd(nightRatio(~isinf(nightRatio)))./meanNightRatio(mySP,1);

            pVocalVec = (dayLightPercBin.*meanDayRatio(mySP,1))+((1-dayLightPercBin).*meanNightRatio(mySP,1));
            pVocalVecCV = (dayLightPercBin*cvDayRatio(mySP,1))+((1-dayLightPercBin)*cvNightRatio(mySP,1));


            dayBins(mySP) = sum(dayPosSet);
            nightBins(mySP) = sum(nightPosSet);
            %NOTE: this pVocal is only used for delphinds. It is hard-coded
            %in density_calc for other species

            inWindow = find(myWeekBinsAllSort>=overallStart&myWeekBinsAllSort<=overallEnd);
            [densityEst,densityCV,siteOrder, spOrder,fnRate,fpRate,pDet, ...
                pVocalTemp,thisGroupSize] = density_calc_report(weekMean(inWindow),weekEffort(inWindow),...
                myWeekBinsAllSort(inWindow),mySpID{mySP},siteCode,[pVocalVec',pVocalVecCV']);

            pVocal(mySP,:) = pVocalTemp; % if used something different internally, adjust here.
            densityTimes = myWeekBinsAllSort(inWindow)';
            if sum(isnan(densityEst))~=length(densityEst)

                [annualChange,annualChangeC1,annualChangeC2,deseasTrendPlusMean,...
                    startValMean,startValCV,endValMean,endValCV,bootSlopeMeanCV,percentChange_10yr]=...
                    Deseason_data_report([densityTimes,densityEst,densityCV],...
                    outDir,[mySpID{mySP},'_',siteCode],[],0,siteCode,[]);
                nktTkt(mySP,1) = 100*(sum(myTSbinAllSort(:,mySP)>=minClickVal)./length(myTimeBinsAllSort));
                nkt(mySP,1) = sum(myTSbinAllSort(:,mySP));

                annualChangeMat(mySP,:) = [nanmean(densityEst),deseasTrendPlusMean(1),...
                    annualChange,annualChangeC1,annualChangeC2];
                stdDensity(mySP,1) = nanstd(densityEst);
                meanDensity(mySP,1) = nanmean(densityEst);
                densityTimesAll{mySP,1} = densityTimes;
                densityEstAll{mySP,1} = densityEst;
                densityCVAll{mySP,1} = densityCV;
                deseasTrendPlusMeanAll{mySP,1} = deseasTrendPlusMean;
                fpRateAll(mySP,1) = fpRate;
                fnRateAll(mySP,1) = fnRate;
                pDetAll(mySP,:) = pDet;
                pVocalAll(mySP,:) = pVocal(mySP,:);
                grpSizeAll(mySP,:) = thisGroupSize;
                startValAll(mySP,:) = [startValMean,startValCV];
                endValAll(mySP,:) = [endValMean,endValCV];
                bootSlopeAll(mySP,:) = bootSlopeMeanCV;
                percentChange_10yrAll(mySP,:) = percentChange_10yr;
            else
                densityEstAll{mySP,1} = {};
                densityCVAll{mySP,1} = {};
            end
        end
        % diel plot
        if (overallEnd-overallStart)>500
            DateTickInterval = round(365/4);
        else
            DateTickInterval = 15;
        end
        if plotsOn
            figure(301);clf
            visPresence([sSet(1:end-1)',sRise(2:end)'],'Color', 'black', ...
                'LineStyle', 'none', 'Transparency', .15, ...
                'Resolution_m', 5/60,'DateTickInterval',DateTickInterval,...
                'DateRange',[overallStart,overallEnd], 'HourTickInterval',...
                2);

            % set anything outside of effort to zero
            myDates = myTimeBinsAllSort(myTSbinAllSort(:,mySP,1)>=minClickVal);
            myCounts = myTSbinAllSort(myTSbinAllSort(:,mySP,1)>=minClickVal,mySP,1);
            keepDates = zeros(size(myDates));
            for iDate = 1:size(effort,1)
                set1 = find(myDates>=effort(iDate,1)&myDates<=effort(iDate,2));
                keepDates(set1)=1;
            end
            myDates(keepDates==0)=[];
            if overallEnd-overallStart<500
                myLineStyle = 'none';
            else
                myLineStyle = '-';
            end

            sightH = visPresence(myDates','Counts',myCounts, ...
                'Resolution_m', 5,'DateTickInterval',DateTickInterval,...
                'LineStyle', 'none', 'DateRange',[overallStart,overallEnd],...
                'Effort', allEffort4Diel, 'Label', outStr,...
                'HourTickInterval', 2,'LineStyle',myLineStyle);
            set(gca, 'YDir', 'reverse','fontweight','bold');  %upside down plot


            box on
            t = annotation('textbox');
            t.String = siteCode;
            if iSite>8
                set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[1 1 2.25 8])
                t.Position = [.65,.83,.05,.1];
            else
                set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[1 1 2.7 8])
                t.Position = [.73,.83,.05,.1];
            end
            t.LineStyle = 'none';
            t.FontSize = 12; t.FontWeight = 'bold';

            %%

            print(gcf,'-dpng','-r300',fullfile(outDir, [mySpID{mySP},'_',siteCode,'_noNorm_weekly_bin_diel.png']))
            %
            saveas(gcf,fullfile(outDir, [mySpID{mySP},'_',siteCode,'_noNorm_weekly_bin_diel.fig']))

        end




    end
    annualChangeMat(:,6) = -annualChangeMat(:,2)./annualChangeMat(:,3);
    %     siteTable = array2table(annualChangeMat,'VariableNames',...
    %         cellstr({'Density';'Yintercept';'Slope';'perc25';'perc75';'TTG'}),'RowNames',mySpID);
    save(fullfile(outDir,[siteCode,'_DensityStats.mat']),'densityTimesAll',...
        'densityEstAll','densityCVAll','mySpID','annualChangeMat','meanDensity','deseasTrendPlusMeanAll',...
        'dayBins','nightBins','pVocal','meanDayRatio','cvDayRatio','allTScountsEffortAdjust','allEffort',...
        'stdDensity','nktTkt','spOrder','fnRateAll','fpRateAll','pDetAll','pVocalAll','myTimeBinsAllSort',...
        'grpSizeAll','startValAll','endValAll','mergedPosDays','mergedDays', 'mergedSpPosDays','percentChange_10yrAll',...
        'bootSlopeAll','-v7.3');

end



