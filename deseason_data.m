function [annualChange,annualChangeC1,annualChangeC2,deseasTrendPlusMean,...
    startValMean,startValCV,endValMean,endValCV,bootSlopePrctile,percentChange_10yr]=...
    Deseason_data_report(myData,savePath,outputFileName,yMax,oilFlag, siteCode,yMax91)
% deseason longterm data
% https://www.mathworks.com/help/econ/parametric-trend-estimation.html
% JAH 12/30/2016

% myData = [matlabDates,meanMeasurement,standard deviation (optional)]% a column matrix
% example: myData = [734139, 1.0, 0.6;
%                    734140, 0.7, 0.4]
% meanMeasurement = mean density of animals, or percent of time present,
% etc., that goes with your dates

% savePath = where you want to save figs, string, example: 'D:\myFigs'
% outputFileName = string to use for naming figs, example: 'MC_Rissos_'

% kef 9/25/2018

% clearvars
monthlyBins = 0; % set to 1 if you want to convert to monthly bins.
% If not, we assume you are providing weekly bins.
savePlots = 1;% Set to 1 if you want to save plots
%% 

outputFileName = [outputFileName];
%% Begin calculations

if size(myData,2)<2
    disp('ERROR: Expecting data at least 2 columns. Check data and delimiter.')
    return
else
    dateVec = myData(:,1);
    meanVec = myData(:,2);
    if size(myData,2)>2
        cvVec = myData(:,3);
    else
        cvVec = nan(size(dateVec));
    end
end

% Test to make sure that the user has supplied weekly bins if they want to
% use weekly bins for estimating seasonality. If the test fails, use monthly
% bins, because that option will re-bin the data anyway.
% if mean(diff(dateVec))~= 7 && monthlyBins == 0
%     fprintf(['WARNING: your data does not appear to be in weekly bins\n',...
%         'Forcing use of monthly bins to handle unknown interval.\n'])
%     monthlyBins = 1;
% end

%% Identify placeholders (NaN or -1), and replace with NaN
missingDataRows = union(find(meanVec<0),find(isnan(meanVec)));
goodRows = setdiff(1:length(meanVec),missingDataRows);

% dateVec(missingDataRows) = NaN;
meanVec(missingDataRows) = NaN;
cvVec(missingDataRows) = NaN;

%% Plot initial data
% figure(91);clf
% plot(dateVec, meanVec,'.');
% set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[1 2 10 5])
% yMax91 = get(gca,'yLim');
% repY91 = repmat(yMax91(2),size(dateVec,1),1);
% bar(dateVec(isnan(meanVec)),repY91(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
%     'EdgeColor',[0.8,0.8,0.8])
% hold on
% plot(dateVec, meanVec,'.');
% xlim([min(dateVec),max(dateVec)])
% datetick('x','mmm ''yy','keepLimits')
% xlabel('Date (month, year)','FontSize',12)
% ylabel('Mean','FontSize',12)
% simpleFName = strrep(strrep(outputFileName,'_','\_'),'.txt','');
% title(simpleFName,'FontSize',10)
% ylim([0,max(repY91)])

%% Call Theil-Sen to estimate slope of entire dataset
dataMat = [dateVec,meanVec];
[estSlope,estOffSet] = TheilSen(dataMat);
trend = estSlope * dateVec + estOffSet;
meanMinusTrend = meanVec - trend;

%% Add trend line to original plot
% figure(91); hold on
% plot(dateVec,trend,'-r')
% hold off
% set(gca,'layer','top')
% legend({'Recording Gaps','Original Data', 'Trend 1: incl. season'})
% if savePlots
%     disp('Saving plots')
%     figName91 = fullfile(savePath,[outputFileName,'_orig_timeseries']);
%     saveas(91,figName91,'fig')
%     print(91,'-dpng','-r600',[figName91,'.png'])
% else
%     disp('No plots saved.')
% end
%% Make second plot with detrended timeseries
% figure(92);clf
% plot(dateVec,meanMinusTrend,'.')
% set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[2 2 10 5])
% datetick('x','mmm ''yy','keepLimits')
% xlabel('Date (month, year)','FontSize',12)
% ylabel('Mean','FontSize',12)
% title(['Detrended Timeseries: ', simpleFName],'FontSize',10)
% xlim([min(dateVec),max(dateVec)])

%% Estimate seasonal component
% How many weeks of data?
%if ~monthlyBins
    fprintf('Calculating seasonal trend using input bins (assumes weekly)\n')
    dateList = weeknum(dateVec);
    onesAndZeros = dummyvar(dateList);
    % remove NaN rows. Could probably do this implicitly.
    onesAndZeros(missingDataRows,:) = [];
    meanMinusTrend(missingDataRows,:) = [];
    seasonalAdjust = onesAndZeros\meanMinusTrend;
    seasonalAdjust(isinf(seasonalAdjust))=0;
    seasonalComponent = onesAndZeros*seasonalAdjust;
    
%else
    fprintf('Calculating seasonal trend using MONTHLY bins\n')
    % this date math is ugly because we are trying to find all of the first
    % of the months between the start and end of the data, including the
    % month that starts before the first point. Since months have
    % different numbers of days, it gets messy. Using financial toolbox
    % functions.
    m_startDate = datenum(datevec(eomdate(min(dateVec)))-[0,1,0,0,0,0])-1;
    m_endDate = eomdate(max(dateVec));
    m_monthStarts = unique(eomdate(m_startDate:15:m_endDate)+1)';
    % datestr(monthStarts) % <-sanity check: now we have datenums for
    % firsts of the months, including start month.
    m_dateList = month(m_monthStarts); % function for figuring out month
    
    % bin the input data using the monthly starts
    [m_binCounts, m_binIdx] = histc(dateVec,m_monthStarts);
    % compute mean for each month
    m_monthlyMeanVec = nan(size(m_monthStarts));
    for iB = 1:length(m_binCounts)
        m_monthlyMeanVec(iB) = nanmean(meanVec(m_binIdx==iB));
    end
    
    % calculate monthly version of trend
    m_monthTrend = estSlope * m_monthStarts + estOffSet;
    m_meanMinusTrend = m_monthlyMeanVec - m_monthTrend;
    m_goodRowsMonthly = find(~isnan(m_monthlyMeanVec));
    m_missingDataRows = find(isnan(m_monthlyMeanVec));
    m_onesAndZeros = dummyvar(m_dateList);
    % remove NaN rows. Could probably do this implicitly.
    m_onesAndZeros(m_missingDataRows,:) = [];
    m_seasonalAdjust = m_onesAndZeros\m_meanMinusTrend(m_goodRowsMonthly);
    m_seasonalComponent = m_onesAndZeros*m_seasonalAdjust;
    
%end

%% Plot seasonal adjustment
figure(93);clf
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 2 5 5])
stem(1:12,m_seasonalAdjust,'ob')


%if monthlyBins

    xlabel('Months','FontSize',12)
    set(gca,'xtick',1:12,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sept','Oct','Nov','Dec'})
    xlim([0.5,12.5])
% else
%     xlabel('Weeks','FontSize',12)
%     xlim([0.5,52.5])
% end
ylabel('Seasonal Adjustment','FontSize',12)
%title(['Estimated Seasonality: ', simpleFName],'FontSize',10)
t = annotation('textbox');
t.String = siteCode;
t.Position = [.82,.8,.05,.1];t.LineStyle = 'none';
t.FontSize = 12; t.FontWeight = 'bold';
set(gca,'fontweight','bold')
if savePlots
    figName93 = fullfile(savePath,[outputFileName,'_seasonality']);
    saveas(93,figName93,'fig')
    print(93,'-dpng','-r300',[figName93,'.png'])
end
% 
% figure(930);clf
% barVals = seasonalAdjust./nanmean(meanVec(goodRows));
% set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[3 2 5 5])
% ax930 = gca;
% % 
% hb = bar(barVals,'Facecolor',[.5,.5,.5]);
% hold on
% yMinVal = ax930.YLim(1);
% yMaxVal = ax930.YLim(2);
% 
% if oilFlag
%     if strcmp(siteCode,'MC')
%         r = rectangle('Position',[4.7143,yMinVal,(12*(90/365)),abs(yMinVal)+yMaxVal],'FaceColor',...
%             [1, 0.6980, 0.6980],'EdgeColor',[1, 0.6980, 0.6980]);
%     elseif strcmp(siteCode,'MCFall')
%         r = rectangle('Position',[9,yMinVal,(12*(90/365)),abs(yMinVal)+yMaxVal],'FaceColor',...
%             [1, 0.6980, 0.6980],'EdgeColor',[1, 0.6980, 0.6980]);
%     elseif strcmp(siteCode,'GCscene')
%         r = rectangle('Position',[4.7143,yMinVal,(12*(90/365)),abs(yMinVal)+yMaxVal],'FaceColor',...
%             [1, 0.6980, 0.6980],'EdgeColor',[1, 0.6980, 0.6980]);
%     elseif strcmp(siteCode,'DTscene')
%         r = rectangle('Position',[4.7143,yMinVal,(12*(90/365)),abs(yMinVal)+yMaxVal],'FaceColor',...
%             [1, 0.6980, 0.6980],'EdgeColor',[1, 0.6980, 0.6980]);
%     end
% end
% %bar(barVals,'Facecolor',[.5,.5,.5]);
% box on
% set(gca,'layer','top')
% yExtent = get(gca,'YLim');
% if oilFlag
%     uistack(r,'down')
% end
% ylim([yMinVal,yMaxVal])
% if monthlyBins
%     xlabel('Months','FontSize',12)
%     xlim([0.5,12.5])
% else
%     xlabel('Weeks','FontSize',12)
%     xlim([0.5,52.5])
% end
% ylabel('Seasonal adjustment','FontSize',12)
% % title([simpleFName,'\n'],'FontSize',10)
% set(gca,'fontsize',18)
% %set(gca,'YLim',yExtent);
% 
% ax930.Position = [.18,.19,.76,.76];

% if savePlots
%     figName930 = fullfile(savePath,[outputFileName,'_seasonality_norm']);
%     saveas(930,figName930,'fig')
%     print(930,'-dtiff','-r300',[figName930,'.tif'])
% end
%% Plot seasonal component of timeseries
% figure(94);clf
% if monthlyBins
%     plot(monthStarts(goodRowsMonthly),seasonalComponent,'.')
% else
%     plot(dateVec(goodRows),seasonalComponent,'.')
% end
% set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[4 2 10 5])
% datetick('x','mmm ''yy','keepLimits')
% xlabel('Date (month, year)','FontSize',12)
% ylabel('Mean','FontSize',12)
% title(['Seasonal component estimate: ',simpleFName], 'FontSize',10)
% xlim([min(dateVec),max(dateVec)])

%% Deseason orignal data and re-estimate slope
if monthlyBins
    deseasonMean = monthlyMeanVec(goodRowsMonthly) - seasonalComponent;
    deseasDataMat = [monthStarts(goodRowsMonthly),deseasonMean];
    origDataMat = [monthStarts(goodRowsMonthly),monthlyMeanVec(goodRowsMonthly)];
else
    deseasonMean = meanVec(goodRows) - seasonalComponent;
    deseasDataMat = [dateVec(goodRows),deseasonMean];
    origDataMat = [dateVec(goodRows),meanVec(goodRows)];
end


[deseasEstSlope,deseasEstOffSet,bootSlope,bootOffset,...
    mBootstrap_Mean_SE,bBootstrap_Mean_SE] = TheilSen(deseasDataMat);
deseasTrend = mBootstrap_Mean_SE(1) * dateVec + bBootstrap_Mean_SE(1);
deseasConf1 = (mBootstrap_Mean_SE(1)-mBootstrap_Mean_SE(2)) * dateVec + (bBootstrap_Mean_SE(1));
deseasConf2 = (mBootstrap_Mean_SE(1)+mBootstrap_Mean_SE(2)) * dateVec + (bBootstrap_Mean_SE(1));

[NodeseasEstSlope,NodeseasEstOffSet,NoDeseas_bootSlope,NoDeseas_bootOffset,...
    NoDeseas_mBootstrap_Mean_SE,NoDeseas_bBootstrap_Mean_SE] = TheilSen(origDataMat);
NoDeseasTrend = NoDeseas_mBootstrap_Mean_SE(1) * dateVec + NoDeseas_bBootstrap_Mean_SE(1);
NoDeseasConf1 = (NoDeseas_mBootstrap_Mean_SE(1)-NoDeseas_mBootstrap_Mean_SE(2)) * dateVec + (NoDeseas_bBootstrap_Mean_SE(1));
NoDeseasConf2 = (NoDeseas_mBootstrap_Mean_SE(1)+NoDeseas_mBootstrap_Mean_SE(2)) * dateVec + (NoDeseas_bBootstrap_Mean_SE(1));
% deseas95perc = (confIntSlope'*dateVec') + repmat(confIntOffset,1,size(dateVec,1));
% annualChange = 100*((deseasTrendMinusMean(end)-deseasTrendMinusMean(1))...
%    )./((dateVec(end)-dateVec(1))./365);% mean(monthlyMeanVec(goodRows))
[~,startDataIdx] = nanmin(dateVec);
[~,endDataIdx] = nanmax(dateVec);

if monthlyBins
    
    deseasTrendPlusMean = mean(monthlyMeanVec(goodRowsMonthly) )+deseasTrend-mean(deseasTrend);
    NodeseasTrendPlusMean = mean(monthlyMeanVec(goodRowsMonthly) )+NoDeseasTrend-mean(NoDeseasTrend);
    
    deseasConf1PlusMean = (mean(monthlyMeanVec(goodRowsMonthly) )+deseasConf1-mean(deseasConf1));
    deseasConf2PlusMean = mean(monthlyMeanVec(goodRowsMonthly) )+deseasConf2-mean(deseasConf2);
    meanVal = nanmean(monthlyMeanVec);
else
    deseasTrendPlusMean = mean(meanVec(goodRows) )+deseasTrend-mean(deseasTrend);
    NodeseasTrendPlusMean = mean(meanVec(goodRows) )+NoDeseasTrend-mean(NoDeseasTrend);

    deseasConf1PlusMean = (mean(meanVec(goodRows) )+deseasConf1-mean(deseasConf1));
    deseasConf2PlusMean = mean(meanVec(goodRows) )+deseasConf2-mean(deseasConf2);
    meanVal = nanmean(meanVec);
    
end

annualChange = (mBootstrap_Mean_SE(1))*365;
NoDeseas_annualChange = (NoDeseas_mBootstrap_Mean_SE(1))*365;
% (100*(deseasTrendPlusMean(endDataIdx)-deseasTrendPlusMean(startDataIdx))...
%     ./deseasTrendPlusMean(startDataIdx))...
%     /((dateVec(endDataIdx)-dateVec(startDataIdx))./365);
annualChangeC1 = (mBootstrap_Mean_SE(1)-(100*mBootstrap_Mean_SE(2)))*365;%
% 100*(deseasConf1PlusMean(endDataIdx)-deseasConf1PlusMean(startDataIdx))...
%     ./deseasConf1PlusMean(startDataIdx)...
%     /((dateVec(endDataIdx)-dateVec(startDataIdx))./365);
annualChangeC2 =  (mBootstrap_Mean_SE(1)+(100*mBootstrap_Mean_SE(2)))*365;
% 100*(deseasConf2PlusMean(endDataIdx)-deseasConf2PlusMean(startDataIdx))...
%     ./deseasConf2PlusMean(startDataIdx)...
%     /((dateVec(endDataIdx)-dateVec(startDataIdx))./365);
startDate = datenum([2010 5 1 0 0 0]);
endDate = datenum([2020 5 1 0 0 0]);
startValMean = mean((bootSlope*startDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend));
startValAll =  (bootSlope*startDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend);

startValStd = std((bootSlope*startDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend));
startValCV = startValStd./startValMean;
endValMean = mean((bootSlope*endDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend));
endValAll = (bootSlope*endDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend);

percentChange_10yr = prctile(100*((endValAll-startValAll)./startValAll),[25,50,75]);
endValStd = std((bootSlope*endDate) + bootOffset + mean(meanVec(goodRows))-mean(deseasTrend));
endValCV = endValStd./endValMean;
bootSlopePrctile = prctile(bootSlope*365,[25,50,75]);
NoDeseas_bootSlopePrctile = prctile(NoDeseas_bootSlope*365,[25,50,75]);

%bootSlopeMeanCV = [mean(bootSlope),std(bootSlope)./mean(bootSlope)];
1;
% annualChange = 100*((endVal-meanVal)./meanVal);

%% Plot deseasoned data and slope
figure(95);clf
subplot(4,1,1)
plot(dateVec, meanVec,'.b');
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[1 2 5 5])
yMax91 = get(gca,'yLim');
repY91 = repmat(yMax91(2),size(dateVec,1),1);
bar(dateVec(isnan(meanVec)),repY91(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
hold on
plot(dateVec, meanVec,'.b');
xlim([min(dateVec),max(dateVec)])
datetick('x','mmm ''yy','keepLimits')
%xlabel('Date (month, year)','FontSize',12)
ylabel({'Weelky Density', 'ind./1000 km^2'},'FontSize',10)
simpleFName = strrep(strrep(outputFileName,'_',' @ '),'.txt','');
title(simpleFName,'FontSize',10)
ylim([0,max(repY91)])
grid on;box on
lowBars = ones(size(repY91));

subplot(4,1,2)
bar(dateVec(isnan(meanVec)),1.05*max(seasonalComponent)*lowBars(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
hold on
bar(dateVec(isnan(meanVec)),1.05*min(seasonalComponent)*lowBars(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
if monthlyBins
    plot(monthStarts(goodRowsMonthly),seasonalComponent,'.b')
else
    plot(dateVec(goodRows),seasonalComponent,'.b')
end
hold off
set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[4 2 5 5])
datetick('x','mmm ''yy','keepLimits')
%xlabel('Date (month, year)','FontSize',12)
ylabel({'Seasonal','Adjustment'},'FontSize',10)
%title(['Seasonal component estimate: ',simpleFName], 'FontSize',10)
xlim([min(dateVec),max(dateVec)])
ylim([min(min(seasonalComponent),0)*1.05,max(max(seasonalComponent),1)*1.05])
box on;grid on

subplot(4,1,3)

bar(dateVec(isnan(meanVec)),1.05*max(deseasonMean)*lowBars(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
hold on
bar(dateVec(isnan(meanVec)),1.05*min(deseasonMean)*lowBars(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
% plot(dateVec,deseas95perc','--k')
if monthlyBins
    plot(monthStarts(goodRowsMonthly),deseasonMean,'.b')
else
    plot(dateVec(goodRows),deseasonMean,'.b')
end

plot(dateVec,deseasTrend','-r')
hold off
%set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[5 2 10 5])
datetick('x','mmm ''yy','keepLimits')
%xlabel('Date (month, year)','FontSize',12)
ylabel({'Deseasoned','Weekly Density'},'FontSize',10)
% title(sprintf('Deseasoned data and slope estimate: %s\n Estimated rate of change: %.3f animals/area/year',...
%     simpleFName,annualChange), 'FontSize',10) 
%legend({'Deseasoned Data', 'Trend 2: de-seasoned'})
xlim([min(dateVec),max(dateVec)])
ylim([min(min(deseasonMean),0)*1.05,max(max(deseasonMean),1)*1.05])
box on;grid on

t3 = annotation('textbox');
t3.String = sprintf('m = %0.2f ind./1000 km^2/year',annualChange);
t3.Position = [.6,.46,.4,.03];t3.LineStyle = 'none';
t3.FontSize = 8; t3.FontWeight = 'bold';


% if savePlots


subplot(4,1,4)
bar(dateVec(isnan(meanVec)),repY91(isnan(meanVec)),1,'FaceColor',[0.8,0.8,0.8],...
    'EdgeColor',[0.8,0.8,0.8])
hold on
if monthlyBins
    plot(monthStarts(goodRowsMonthly),monthlyMeanVec(goodRowsMonthly),'.')
else
    plot(dateVec(goodRows),meanVec(goodRows),'.b')
end
plot(dateVec,deseasTrendPlusMean','-r')
plot(dateVec,deseasConf1PlusMean','--r')
plot(dateVec,deseasConf2PlusMean','--r')
hold off
% plot(dateVec,deseas95perc','--k')
%set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[5 2 10 5])
datetick('x','mmm ''yy','keepLimits')
xlabel('Date (month, year)','FontSize',10)
ylabel({'Weekly Density',' ind./1000 km^2'},'FontSize',10)
% title(sprintf('Deseasoned data and slope estimate: %s\n Estimated rate of change: %.3f animals/area/year',...
%     simpleFName,annualChange), 'FontSize',10)
%legend({'Deseasoned Data', 'Trend 2: de-seasoned'})
xlim([min(dateVec),max(dateVec)])
box on;grid on
% t4 = annotation('textbox');
% t4.String = sprintf('m = %0.2f [%0.2f,%0.2f]', NoDeseas_bootSlopePrctile(1),...
%     NoDeseas_bootSlopePrctile(2),NoDeseas_bootSlopePrctile(3));
% t4.Position = [.60,.24,.3,.03];t4.LineStyle = 'none';
% t4.FontSize = 8; t4.FontWeight = 'bold';%
    figName95 = fullfile(savePath,[outputFileName,'_deseasoned']);
    saveas(95,figName95,'fig')
    print(95,'-dpng','-r300',[figName95,'.png'])

%% Plot original data and slope
figure(96);clf
if isempty(yMax)
    cvMax = max(nanmax(cvVec),0);
else
    cvMax = yMax;
end

if ~isnan(nanmax(cvVec)) % have cvs
    hE = errorbar(dateVec, meanVec,cvVec,'.k','MarkerSize',1,'MarkerFaceColor',[.5 .5 .5]);
    hE.CapSize = 3;
    
    hold on
    if (nanmax(dateVec)-nanmean(dateVec))>500
        myEdgeColor = 'b';
    else
        myEdgeColor = 'k';
    end
    bar(dateVec, meanVec,1,'FaceColor',[.5 .5 .5],'EdgeColor',myEdgeColor)
    if isempty(yMax)
        yMax = nanmax(meanVec+cvVec)*1.1;
    end
else
    %no cvs
    plot(dateVec,meanVec,'ok','MarkerSize',4,'MarkerFaceColor',[.5 .5 .5]);
    if isempty(yMax)
        yMax = max(nanmax(meanVec)*1.1,1);
    end
end
xlim([min(dateVec),max(dateVec)])
ylim([0,max(yMax,1)])
set(gca,'FontSize',10)
hold on
% plot(dateVec,deseasTrendPlusMean,'--r','lineWidth',2)
% plot(dateVec,deseasConf1PlusMean,'--r','lineWidth',2)
% plot(dateVec,deseasConf2PlusMean,'--r','lineWidth',2)
grid on
if isempty(yMax91)    
    yMax91 = get(gca,'yLim');
end
repY91 = repmat(yMax91(2),size(dateVec,1),1);
bar(dateVec(isnan(meanVec)),1.1*(cvMax+repY91(isnan(meanVec))),...
    1,'FaceColor',[0.76,0.74,0.99],...
    'EdgeColor',[0.76,0.74,0.99])
hold on
set(gcf,'units','inches','PaperPositionMode','auto')
if (max(dateVec)-min(dateVec))<500
    set(gcf,'OuterPosition',[6 2 5 5])
else
    set(gcf,'OuterPosition',[6 2 5 3.5])
end
datetick('x','mmm ''yy','keepLimits')
xlabel('Date (month/year)','FontSize',10)
ylabel({'Weekly Mean Density';'(animals/1000 km^2)'},'FontSize',1)
% legend({'Recording Gaps','Original Data', 'Trend 2: de-seasoned'})
%title(siteCode,'FontSize',10)
xlim([min(dateVec),max(dateVec)])
tickList = get(gca,'xtick');
thisYear = year(tickList)';
thisMonth = month(tickList)';
thisDay = ones(size(thisMonth));
thisStartMonth = datenum([thisYear,thisMonth,thisDay]);
set(gca,'xtick',thisStartMonth,'fontsize',12,'fontweight','bold')
datetick('x','mm/yy','keeplimits');grid on

set(gca,'layer','top')
if oilFlag
    plot([datenum([2010,4,20]),datenum([2010,7,19])],[yMax,yMax]*.95,...
        'Color',[.64,.08,.18],'LineWidth',4)
end
t = annotation('textbox');
t.String = siteCode;
t.Position = [.8,.8,.05,.1];t.LineStyle = 'none';
t.FontSize = 12; t.FontWeight = 'bold';
if savePlots
    figName96 = fullfile(savePath,[outputFileName,'_longTermTrend','_300dpi']);
    saveas(96,figName96,'fig')
    print(96,'-dpng','-r300',[figName96,'.png'])
end

% confIntSlopeAnnual = confIntSlope*(365)./mean(diff(dateVec));
% endValConf1 = meanVal+(confIntSlope(1)*(365));
% endValConf2 = meanVal+(confIntSlope(2)*(365));
% annualChangeConf1 = 100*((endValConf1-meanVal)./meanVal);
% annualChangeConf2 = 100*((endValConf2-meanVal)./meanVal);
%%
figure(98);clf
if ~isnan(nanmax(cvVec)) % have cvs
    hE = errorbar(dateVec, meanVec,cvVec,'.k','MarkerSize',1,'MarkerFaceColor',[.5 .5 .5]);
    hE.CapSize = 1;
    
    hold on
    if (nanmax(dateVec)-nanmean(dateVec))>500
        myEdgeColor = 'b';
    else
        myEdgeColor = 'k';
    end
    %bar(dateVec, meanVec,1,'FaceColor',[.5 .5 .5],'EdgeColor',myEdgeColor)
    plot(dateVec,meanVec,'ok','MarkerSize',4,'MarkerFaceColor',[.5 .5 .5]);

    if isempty(yMax)
        yMax = nanmax(meanVec+cvVec)*1.1;
    end
else
    %no cvs
    plot(dateVec,meanVec,'ok','MarkerSize',4,'MarkerFaceColor',[.5 .5 .5]);
    if isempty(yMax)
        yMax = nanmax(meanVec)*1.1;
    end
end
xlim([min(dateVec),max(dateVec)])
try ylim([0,yMax]) 
end
set(gca,'FontSize',11)
hold on 
%plot(dateVec,deseasTrendPlusMean,'--r','lineWidth',2)
%plot(dateVec,deseasConf1PlusMean,'--r','lineWidth',2)
%plot(dateVec,deseasConf2PlusMean,'--r','lineWidth',2)
% title(sprintf('%.2f animals/1000km^2 per year',...
%     annualChange),'FontSize',10)
grid on
if isempty(yMax91)
    yMax91 = get(gca,'yLim');
end
repY91 = repmat(yMax91(2),size(dateVec,1),1);
bar(dateVec(isnan(meanVec)),1.1*(cvMax+repY91(isnan(meanVec))),...
    1,'FaceColor',[0.76,0.74,0.99],...
    'EdgeColor',[0.76,0.74,0.99])
hold on
%set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[6 2 10 4])
if (max(dateVec)-min(dateVec))<500
    set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[6 2 5 5])
else
    set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[6 2 5 3.5])
end
datetick('x','mmm ''yy','keepLimits')
xlabel('Date (month/year)','FontSize',12)
ylabel({'Weekly Mean Density';'(animals/1000 km^2)'},'FontSize',12)
% legend({'Recording Gaps','Original Data', 'Trend 2: de-seasoned'})
%title(siteCode,'FontSize',10)
xlim([min(dateVec),max(dateVec)])
tickList = get(gca,'xtick');
thisYear = year(tickList)';
thisMonth = month(tickList)';
thisDay = ones(size(thisMonth));
thisStartMonth = datenum([thisYear,thisMonth,thisDay]);
set(gca,'xtick',thisStartMonth,'fontsize',10,'fontweight','bold')
datetick('x','mm/yy','keeplimits');grid on

set(gca,'layer','top')
if oilFlag
    plot([datenum([2010,4,20]),datenum([2010,7,19])],[yMax,yMax]*.95,...
        'Color',[.64,.08,.18],'LineWidth',4)
end
t = annotation('textbox');
t.String = siteCode;
t.Position = [.8,.8,.05,.1];t.LineStyle = 'none';
t.FontSize = 12; t.FontWeight = 'bold';
set(98,'PaperPosition',get(96,'PaperPosition'))
if savePlots
    figName98 = fullfile(savePath,[outputFileName,'_longTermTrendwLine','_300dpi']);
    saveas(98,figName98,'fig')
    print(98,'-dpng','-r300',[figName98,'.png'])
end
%%
disp(outputFileName)

fprintf('Estimated annual rate of change is %.3f, [25,75] confidence interval = [%.3f %.3f] \n',...
    annualChange,annualChangeC1,annualChangeC2)


