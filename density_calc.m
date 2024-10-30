function [densityEst,densityCV,siteOrder, spOrder,thisfnRate,thisfpRate, ...
    thispDet,pVocalNew,thisGroupSize] = density_calc_report(binCount,totalBinCount,...
    binTimes,speciesCode,siteCode,pVocal)
% load("E:\Code\Gulf of Mexico\pVocal_GOM.mat")
spOrder = {'Pm','Zc','Me','Md','Gg','UD','UD_LF','Kspp'};
siteOrder = {'AC','CE','DC','DT','GC','LC','MC','MP','MR','Y1B','Y1C','Y1D'};
groupSizeMat = [6.1 0.18; %Pm Richter et al 2008 Social Structure of Sperm whales
    1.95, 0.05;% Zc
    2.26, 0.03;% Me
    1.95, 0.05;% Md
    12.9, 0.46; % Gg
    51.2, 0.30; % UD_HF
    32.8, 0.48; %UD_LF
    1.58,0.09]; %Kspp


if strcmp(speciesCode,'UD_all')
    speciesCode = 'UD';
end

if strcmp(speciesCode,'Pm')  
    pVocal = [.81,0.06];
elseif strcmp(speciesCode,'Zc')
    pVocal = [ 0.254, 0.17];
elseif strcmp(speciesCode,'Me')
    pVocal = [0.471, 0.09];% Me
elseif strcmp(speciesCode,'Md')
    pVocal = [0.165,0.075]; % Warren et al 2017 for Blainvilles
elseif strcmp(speciesCode,'Kspp')
    pVocal = [0.197, 0.121]; %Kspp
end
% elseif strcmp(speciesCode,'Gg')
%     pVocal = [mean([1,myPV(siteCode,:).Gg]),myPVCV(siteCode,:).Gg];
% elseif strcmp(speciesCode,'UD')
%     pVocal = [mean([1,myPV(siteCode,:).UD]),myPVCV(siteCode,:).UD];
% elseif strcmp(speciesCode,'UD_3p')
%     pVocal = [mean([1,myPV(siteCode,:).UD_3p]),myPVCV(siteCode,:).UD_3p];
%     if isnan(pVocal)
%          pVocal = [mean([1,myPV(end,:).UD_3p]),myPVCV(end,:).UD_3p];
%     end
% elseif strcmp(speciesCode,'UD_LF')
%     pVocal = [mean([1,myPV(siteCode,:).UD_LF]),myPVCV(siteCode,:).UD_LF];
%     if isnan(pVocal)
%          pVocal = [mean([1,myPV(end,:).UD_LF]),myPVCV(end,:).UD_LF];
%     end
% end


pVocal(pVocal(:)>1,1)=1;

if size(pVocal,1)==1
    pVocal =repmat(pVocal,length(totalBinCount),1);
end
regIdx = 1:length(binCount);
altIdx = [];

W = [35;4;4;4;5;5;5;1.5];%km

% pDet = [0.592,0.048,0.030,0.012,0.196,0.081,0.081,0.196,0.208;
% 0.522, 0.284, 0.213, 0.005, 0.200, 0.087, 0.087, 0.200, 0.222;
% NaN ,   NaN , NaN  ,   NaN, 0.196, 0.125, 0.125, 0.196, NaN  ;
% 0.387, 0.293, 0.224, 0.235, 0.208, 0.084, 0.084, 0.208, 0.269;% DT
% 0.436, 0.294, 0.225, 0.235, 0.209, 0.084, 0.084, 0.209, 0.267;% GC
% 0.584, 0.294, 0.224, 0.234, 0.190, 0.077, 0.077, 0.190, 0.195;
% 0.377, 0.297, 0.226, 0.237, 0.209, 0.092, 0.092, 0.209, 0.288;% MC
% NaN  , NaN  , NaN  , NaN  , NaN  , 0.1045  , NaN  , NaN  , NaN  ;% MP
% 0.471, 0.302, 0.231, 0.243, 0.203, 0.080, 0.080, 0.203, 0.206;
% 0.567, 0.294, 0.221, 0.224, 0.192, 0.071, 0.071, 0.192, 0.195;
% 0.569, 0.294, 0.223, 0.231, 0.190, 0.077, 0.077, 0.190, 0.195;
% 0.567, 0.294, 0.225, 0.234, 0.207, 0.107, 0.107, 0.207, 0.223];
% 
% 
% pDetAlt = [0.143,0.324,0.248,0.268,0.263,0.120,0.120,0.263,0.242];
% 
% pDetAlt_std = [0.0082,0.027,0.022,0.026,0.0050,0.0033,0.0033,0.0050,0.038];
% pDetAlt_cv = pDetAlt_std./pDetAlt;
% 
% pDet_std = [0.038,0.025,0.019,0.022,0.044,0.024,0.024,0.044,0.021;
% 0.037,0.025,0.021,0.022,0.042,0.028,0.028,0.042,0.021;
% NaN  ,NaN  ,NaN  ,NaN  ,0.024,0.024,0.024,0.024,  NaN;
% 0.033,0.023,0.020,0.023,0.045,0.025,0.025,0.045,0.050;%DT
% 0.027,0.024,0.020,0.022,0.044,0.025,0.025,0.044,0.052;%GC
% 0.052,0.023,0.021,0.022,0.041,0.024,0.024,0.041,0.022;
% 0.021,0.023,0.020,0.022,0.044,0.026,0.026,0.044,0.051;%MC
% NaN  ,NaN  ,NaN  ,NaN  ,NaN  ,0.015,NaN  ,NaN  ,NaN  ;%MP
% 0.027,0.024,0.020,0.022,0.044,0.024,0.024,0.044,0.021;
% 0.049,0.025,0.020,0.022,0.041,0.022,0.022,0.041,0.021;
% 0.053,0.023,0.020,0.022,0.045,0.023,0.023,0.045,0.024;
% 0.049,0.025,0.020,0.022,0.037,0.025,0.025,0.037,0.019];% change to stdError??
% 
%siteOrder = {'AC','CE','DC','DT','GC','LC','MC','MP','MR','Y1B','Y1C','Y1D'};

pDet = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0.2300,0.1353,0.2694,NaN;
    0.2265,0.2926,0.2237,0.2357,0.2111,0.0871,0.3048,0.3142;
    0.205,0.297,0.225,0.235,0.207,0.086,0.301,0.311;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0.1713,0.2971,0.2259,0.2384,0.2125,0.0927,0.2999,0.3333;
    NaN,NaN,NaN,NaN,0.2306,0.1054,0.2984,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];


pDetAlt = [0.09897,0.3222,0.2490,0.2682,0.2698,0.1268,0.3309,0.34499;  %MC_800
    NaN,NaN,NaN,NaN,0.206578770575475,0.0855890007691445,0.260637253365992,NaN%MP_summer
    NaN,NaN,NaN,NaN,0.214063529203445,0.130921970920458,0.252486474912083,NaN;%DC Summer
    0.206481496354950,0.296525322848827,0.225490512374744,0.235443516792955,0.207037073217917,0.0859480473806440,0.300930488938523,0.310993026007223;%GC summer
    0.174248447330193,0.297099831379855,0.225762875370046,0.238410142504594,0.212525207351267,0.0926596898607413,0.299897887348185,0.333314527674074];%MC summer

pDetAlt_std = [0.0149809582622970,0.0266793095698511,0.0216005829540382,0.0259557669945966,0.0538010039898708,0.0297649333873953,0.0688363022001394,0.0272804638002064;  %MC_800   
	NaN,NaN,NaN,NaN,0.0516804545003124,0.0205566298739898,0.0645626341666276,NaN;%MP_summer
    NaN,NaN,NaN,NaN,0.0390059158395684,0.0244355610601596,0.0457212926630235,NaN; %DC summer
    0.0229171214202331,0.0234333253115036,0.0194086175683835,0.0232354565911822,0.0424330167190960,0.0254693705414328,0.0605247642701456,0.0277176691450340;%GC summer
   0.0185871630919222,0.0244968702896185,0.0212546769056669,0.0226579625389126,0.0403761327407216,0.0259485872348779,0.0614036717809850,0.0283542020007612];

pDetAlt_cv = pDetAlt_std./pDetAlt;

pDet_std = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0.0407827914252636,0.0316622675733983,0.0446892707064564,NaN;
    0.0288,0.0247214269940241,0.0196706120331657,0.0223061919243735,0.0438120821185076,0.0249240791696105,0.0646521789359240,0.0287265552599061;
    0.0229,0.0234333253115036,0.0194086175683835,0.0232354565911822,0.0424330167190960,0.0254693705414328,0.0605247642701456,0.0277176691450340;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0.0185,0.024497,0.02125,0.0227,0.0403761327407216,0.0259485872348779,0.0614036717809850,0.0283542020007612;
    NaN,NaN,NaN,NaN,0.0481419416150255,0.0212917703973071,0.0633131814985385,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];


pDet_cv = pDet_std./pDet;
thisSpIdx = find(strcmpi(speciesCode,spOrder));
thisSiteIdx = find(strcmp(siteCode,siteOrder));

errorRate = [4.3, 4.3, 0.1, 2.3, 4.1, 4.3,6.5,NaN,4.3,4.3,4.3,4.3
0.6	0.6	 0.1	1.4	 0.2	0.6	 0.1	NaN	0.6	0.6	0.6	0.6
0.3	0.3	 0.1	0.7	 0.2	0.3	 0.1	NaN	0.3	0.3	0.3	0.3
0.0	0.0	 0.0	0.0	 0.1	0.0  0.0	NaN	0.0	0.0	0.0	0.0
2.5	2.5	 4.6	3.6	 2.0	2.5	 1.9	NaN	2.5	2.5	2.5	2.5
4.9	4.9	14.1	6.1	 4.7	4.9	 3.8	0	4.9	4.9	4.9	4.9
4.9	4.9	14.1	6.1	 4.7	4.9	 3.8	NaN	4.9	4.9	4.9	4.9
0.1	0.1	 0.7	0.1	 0.1	0.1	 0.2	NaN	0.1	0.1	0.1	0.1]/100;
% 
fpRate = [4.0	4.0	0.1	1.0	6.5	4.0	4.4	NaN	4.0	4.0	4.0	4.0
0.2	0.2	0.1	0.5	0.1	0.2	0.1	NaN	0.2	0.2	0.2	0.2
0.1	0.1	0.1	0.2	0.1	0.1	0.1	NaN	0.1	0.1	0.1	0.1
0.0	0.0	0.0	0.0	0.1	0.0	0.0	NaN	0.0	0.0	0.0	0.0
1.9	1.9	4.6	2.5	1.7	1.9	1.4	NaN	1.9	1.9	1.9	1.9
3.3	3.3	8.5	4.6	3.2	3.3	2.1	0	3.3	3.3	3.3	3.3
3.3	3.3	8.5	4.6	3.2	3.3	2.1	NaN	3.3	3.3	3.3	3.3
0.0	0.0	0.7	0.0	0.0	0.0	0.1	NaN	0.0	0.0	0.0	0.0]/100;

fnRate = [5.1 5.1 34.8 5.2 3.0	5.1	7.0	NaN	5.1	5.1	5.1	5.1
7.1  7.1  0.0  7.0	8.3	 7.1  6.1  NaN	7.1	7.1	7.1	7.1
7.4	 7.4  0.0  8.7	5.3	 7.4  8.3  NaN	7.4	7.4	7.4	7.4
33.3 33.3 0.0  100.0 0.0 33.3  0.0 NaN 33.3 33.3 33.3 33.3
15.4 15.4 16.7 7.4	20.6 15.4 18.2 NaN 15.4 15.4 15.4 15.4
11.7 11.7 25.7 10.1 11.3 11.7 13.8 0 11.7 11.7 11.7 11.7
11.7 11.7 25.7 10.1 11.3 11.7 13.8 NaN 11.7 11.7 11.7 11.7
26.5 26.5 66.7 33.3 25.4 26.5 20.7 NaN 26.5 26.5 26.5 26.5]/100;

altIdx = [];
altRow = [];
% if strcmp(siteCode,'MC')&& ~strcmp(speciesCode,'Pm')
%     altDates = [datenum([2014	04	23	00	00	00]),	datenum([2017	05	17	0	00	00])];%mc09-11
% 	altIdx = find(binTimes>=altDates(1)& binTimes<=altDates(2));
%     altRow = 1;
% 
if strcmp(siteCode,'MP')
    altIdx = find(month(binTimes)>=5 & month(binTimes)<=10);% summer months
    altRow = 2;
end

if strcmp(siteCode,'DC')
    altIdx = find(month(binTimes)>=5 & month(binTimes)<=10);% summer months
    altRow = 3;
end

if strcmp(siteCode,'GC')&& max(strcmp(speciesCode,'Gg'),strcmp(speciesCode,'UD_LF'))
    altIdx = find(month(binTimes)>=5 & month(binTimes)<=10);% summer months
    altRow = 4;
end

if strcmp(siteCode,'MC')&& strcmp(speciesCode,'UD_LF')
    altIdx = find(month(binTimes)>=5 & month(binTimes)<=10);% summer months
    altRow = 5;
end
densityNumerator(regIdx,1) = binCount(regIdx).*groupSizeMat(thisSpIdx,1)*(1-fpRate(thisSpIdx,thisSiteIdx))*(1+fnRate(thisSpIdx,thisSiteIdx));
densityNumerator(altIdx,1) = binCount(altIdx).*groupSizeMat(thisSpIdx,1)*(1-fpRate(thisSpIdx,thisSiteIdx))*(1+fnRate(thisSpIdx,thisSiteIdx));
densityDenominator(regIdx,1) = pi*(W(thisSpIdx)^2).*pVocal(regIdx,1).*pDet(thisSiteIdx,thisSpIdx).*(totalBinCount(regIdx));
densityDenominator(altIdx,1) = pi*(W(thisSpIdx)^2).*pVocal(altIdx,1).*pDetAlt(altRow,thisSpIdx).*(totalBinCount(altIdx));

densityEst = (densityNumerator./densityDenominator).*1000;
densityCV(regIdx,1) = sqrt((densityEst(regIdx,1).^2)*((nanmean(pVocal(:,2)).^2)+ ...
    (pDet_cv(thisSiteIdx,thisSpIdx).^2)+ (groupSizeMat(thisSpIdx,2).^2)));
if ~isempty(altIdx)
    densityCV(altIdx,1) = sqrt((densityEst(altIdx,1).^2).*((nanmean(pVocal(:,2)).^2)+ ...
        (pDetAlt_cv(altRow,thisSpIdx).^2)+ (groupSizeMat(thisSpIdx,2).^2)));
end

pVocalNew(1) = nanmean(pVocal(:,1));
pVocalNew(2) = nanmean(pVocal(:,2));

thisfnRate = fnRate(thisSpIdx,thisSiteIdx);% this order of indices is right, backward from the others.
thisfpRate = fpRate(thisSpIdx,thisSiteIdx);% this order of indices is right, backward from the others.
thispDet = [pDet(thisSiteIdx,thisSpIdx),pDet_cv(thisSiteIdx,thisSpIdx)];
thisGroupSize = groupSizeMat(thisSpIdx,:);

