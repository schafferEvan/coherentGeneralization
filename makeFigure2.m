
% RELEVANT MAT FILES WERE GENERATED BY:                      
% calculateSumFigParts.m             (for examples and class odors)
% calculate4thOrderPopResp_v1.m      (for scaling plots, N4=1)
% _flyModel/calculate4thOrderPopResp_Fly_antiHebbSmall.m
% reorganized and moved pop to sup fig on 8/3/17

global classLevel1 classLevel2 greenColorAdjust class1Color class2Color

classLevel1 = 0.7;
classLevel2 = 0.3;
greenColorAdjust = 0.15;
class1Color = [0.0    0.7    0.2]; %[0    0.9    0.2]; %[0, classLevel1+greenColorAdjust, greenColorAdjust];
class2Color = [0    0.5    0.6]; %[0, classLevel2+greenColorAdjust, classLevel2+greenColorAdjust];


load('_simResults/partsForSumFig.mat');
%initializeGlobals;
makeMixOfExpertsFig = false;

%% 2nd Model Figure (fig 3)


paperFig1 = figure; 
set(gcf,'color','w')
set(gcf,'Units','Centimeters')
pos(3)=11.4; %17.4; %8.5;%pos(3)*1.5;
pos(4)=14.0;%pos(4)*1.5;
set(gcf,'Position',pos)

aiAxes = axes('position',[.11  .50  0.2*pos(4)/pos(3)  0.2]); hold on
aiiAxes = axes('position', [.11  .72  0.2*pos(4)/pos(3)  0.2]); hold on
annotation('textbox', [.02 .94 .05 .05],...
    'String', 'A','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');

bAxes = axes('position',[.55  .5  .42  .42]);
annotation('textbox', [.43 .94 .05 .05],...
    'String', 'B','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');

cAxes = axes('position',[.11  .08  .35  .26]); 
annotation('textbox', [.02 .37 .05 .05],...
    'String', 'C','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');

diAxes = axes('position',[.58  .08  .08*pos(4)/pos(3)  .08]);
diiAxes = axes('position',[.58  .26  .08*pos(4)/pos(3)  .08]); 
annotation('textbox', [.49 .37 .05 .05],...
    'String', 'D','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');

annotation('textbox', [.705 .37 .05 .05],...
    'String', 'E','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');


%% mouse 1 vs mouse 2 training -------------  v1

popDir = [pwd,'/_simResults/_pirSims/_bigParts/'];
popDL = dir([popDir,'pop4thOrderResp*']);
popTot = load([popDir,popDL(1).name]);
for j=2:length(popDL)
    popTot(j) = load([popDir,popDL(j).name]);
    popTot(j).rAgr(popTot(j).rAgr==0)=nan; 
    popTot(j).rAgr90(popTot(j).rAgr90==0)=nan; 
    popTot(j).rAcc(popTot(j).rAcc==0)=nan; 
    popTot(j).rCC(popTot(j).rCC==0)=nan; 
    popTot(j).rAgrSingle(popTot(j).rAgrSingle==0)=nan; 
    popTot(j).rAgrSum(popTot(j).rAgrSum==0)=nan; 
    popTot(j).snr(popTot(j).snr==0)=nan; 
    popTot(j).ZCat(popTot(j).ZCat==0)=nan; 
    popTot(j).Z2Cat(popTot(j).Z2Cat==0)=nan;
end


normQt1 = trainedFourthOrder90pClass1(trainedOdorNum);
normQt2 = trainedFourthOrder90pClass1_mouse2(trainedOdorNum);

mouse1class1 = trainedFourthOrder90pClass1/trainedFourthOrder90pClass1(trainedOdorNum); %trainedFourthOrder90pMix;
mouse2class1 = trainedFourthOrder90pClass1_mouse2/trainedFourthOrder90pClass1_mouse2(trainedOdorNum); %trainedFourthOrder90pMix2_mouse2
%mouse1class2 = gCCn(2:length(gCCn),trainedOdorNum)'/max(gCCn(2:length(gCCn),trainedOdorNum)); %trainedFourthOrder90pMix;
mouse1class2 = trainedFourthOrder90pClass2/trainedFourthOrder90pClass1(trainedOdorNum); %trainedFourthOrder90pMix;
mouse2class2 = trainedFourthOrder90pClass2_mouse2/trainedFourthOrder90pClass1_mouse2(trainedOdorNum); %trainedFourthOrder90pMix2_mouse2
mouse1rand = trainedFourthOrder90pRand/trainedFourthOrder90pClass1(trainedOdorNum); %trainedFourthOrder90pMix;
mouse2rand = trainedFourthOrder90pRand_mouse2/trainedFourthOrder90pClass1_mouse2(trainedOdorNum); %trainedFourthOrder90pMix2_mouse2


EFmarkerSize = 9;
axes(aiAxes)
xMax = full(max( [ max(mouse1class1),max(mouse1class2),max(mouse1rand) ] ));
yMax = full(max( [ max(mouse2class1),max(mouse2class2),max(mouse2rand) ] ));
xMin = full(min( [ min(mouse1class1),min(mouse1class2),min(mouse1rand) ] ));
yMin = full(min( [ min(mouse2class1),min(mouse2class2),min(mouse2rand) ] ));

mTh1 = full(median(mouse1class1));
mTh2 = full(median(mouse2class1));
rectangle('Position',[mTh1,yMin,xMax-mTh1,mTh2-yMin],'facecolor',.9*[1,1,1],'edgecolor','none')
rectangle('Position',[xMin,mTh2,mTh1-xMin,yMax-mTh2],'facecolor',.9*[1,1,1],'edgecolor','none')

plot(mouse1rand,mouse2rand,'.','color',[0,0,0],'markersize',EFmarkerSize); %[.7 .7 .7]
plot(mouse1class2,mouse2class2,'.','color',class2Color,'markersize',EFmarkerSize); %[.7 .7 .7]
plot(mouse1class1,mouse2class1,'.','color',class1Color,'markersize',EFmarkerSize); %[.7 .7 .7]

box off
ylim([yMin yMax]);
xlim([xMin xMax]);
fracAgree70 = ( sum( (mouse2class1>mTh2) & (mouse1class1>mTh1)  ) ...
    +  sum( (mouse2class1<mTh2) & (mouse1class1<mTh1)  ) )/length(mouse1class1);
Ath70 = 1/(1-.5)*(fracAgree70-0.5);
mTh1_30 = full(median(mouse1class2));
mTh2_30 = full(median(mouse2class2));
fracAgree30 = ( sum( (mouse2class2>mTh2_30) & (mouse1class2>mTh1_30)  ) ...
    +  sum( (mouse2class2<mTh2_30) & (mouse1class2<mTh1_30)  ) )/length(mouse1class2);
Ath30 = 1/(1-.5)*(fracAgree30-0.5);
mTh1_0 = full(median(mouse1rand));
mTh2_0 = full(median(mouse2rand));
fracAgree0 = ( sum( (mouse2rand>mTh2_0) & (mouse1rand>mTh1_0)  ) ...
    +  sum( (mouse2rand<mTh2_0) & (mouse1rand<mTh1_0)  ) )/length(mouse1rand);
Ath0 = 1/(1-.5)*(fracAgree0-0.5);plot([mTh1,mTh1],ylim,'k','linewidth',1)
plot(xlim,[mTh2,mTh2],'k','linewidth',1) 
% set(gca,'Xtick',[mTh1,xyMax])
% set(gca,'Ytick',[mTh2,xyMax])
set(gca,'Xtick',full([mTh1,xMax]))
set(gca,'Ytick',full([mTh2,yMax]))
set(gca,'XtickLabel',{'\theta','1'});%,'fontsize',11)
set(gca,'YtickLabel',{'\theta','1'});%,'fontsize',11)
xlabel('readout 1','fontsize',12)
ylabel('readout 2','fontsize',12)

display(['trained A_{0.5}=',num2str(.01*round(100*Ath70)),', ',num2str(.01*round(100*Ath30)),', ',num2str(.01*round(100*Ath0))]);
c1 = corrcoef(mouse2class1,mouse1class1);
c2 = corrcoef(mouse2class2,mouse1class2);
c0 = corrcoef(mouse2rand,mouse1rand);
display(['trained CC=',num2str(c1(2)),', ',num2str(c2(2)),', ',num2str(c0(2))]);



%% mouse 1 vs mouse 2 training -------------  v4

mouse1class1 = UNTRAINEDFourthOrder90pClass1/normQt1; %UNTRAINEDFourthOrder90pClass1(trainedOdorNum); %trainedFourthOrder90pMix;
mouse2class1 = UNTRAINEDFourthOrder90pClass1_mouse2/normQt2; %UNTRAINEDFourthOrder90pClass1_mouse2(trainedOdorNum); %trainedFourthOrder90pMix2_mouse2
%mouse1class2 = gCCn(2:length(gCCn),trainedOdorNum)'/max(gCCn(2:length(gCCn),trainedOdorNum)); %trainedFourthOrder90pMix;
mouse1class2 = UNTRAINEDFourthOrder90pClass2/normQt1; %max(UNTRAINEDFourthOrder90pClass1); %trainedFourthOrder90pMix;
mouse2class2 = UNTRAINEDFourthOrder90pClass2_mouse2/normQt2; %max(UNTRAINEDFourthOrder90pClass1_mouse2(2:length(gCCn))); %trainedFourthOrder90pMix2_mouse2
mouse1rand = UNTRAINEDFourthOrder90pRand/normQt1; %UNTRAINEDFourthOrder90pClass1(trainedOdorNum); %trainedFourthOrder90pMix;
mouse2rand = UNTRAINEDFourthOrder90pRand_mouse2/normQt2; %UNTRAINEDFourthOrder90pClass1_mouse2(trainedOdorNum); %trainedFourthOrder90pMix2_mouse2


axes(aiiAxes)
xMax = full(max( [ max(mouse1class1),max(mouse1class2),max(mouse1rand) ] ));
yMax = full(max( [ max(mouse2class1),max(mouse2class2),max(mouse2rand) ] ));
xMin = full(min( [ min(mouse1class1),min(mouse1class2),min(mouse1rand) ] ));
yMin = full(min( [ min(mouse2class1),min(mouse2class2),min(mouse2rand) ] ));

%xyMax = full(tempsort( round(.01*length(tempsort)) )); %this includes 99% of points
mTh1 = full(median(mouse1class1));
mTh2 = full(median(mouse2class1));
rectangle('Position',[mTh1,yMin,xMax-mTh1,mTh2-yMin],'facecolor',.9*[1,1,1],'edgecolor','none')
rectangle('Position',[xMin,mTh2,mTh1-xMin,yMax-mTh2],'facecolor',.9*[1,1,1],'edgecolor','none')

plot(mouse1rand,mouse2rand,'.','color',[0,0,0],'markersize',EFmarkerSize); %[.7 .7 .7]
plot(mouse1class2,mouse2class2,'.','color',class2Color,'markersize',EFmarkerSize); %[.7 .7 .7]
plot(mouse1class1,mouse2class1,'.','color',class1Color,'markersize',EFmarkerSize); %[.7 .7 .7]

box off
ylim([yMin yMax]);
xlim([xMin xMax]);
fracAgree70 = ( sum( (mouse2class1>mTh2) & (mouse1class1>mTh1)  ) ...
    +  sum( (mouse2class1<mTh2) & (mouse1class1<mTh1)  ) )/length(mouse1class1);
Ath70 = 1/(1-.5)*(fracAgree70-0.5);
mTh1_30 = full(median(mouse1class2));
mTh2_30 = full(median(mouse2class2));
fracAgree30 = ( sum( (mouse2class2>mTh2_30) & (mouse1class2>mTh1_30)  ) ...
    +  sum( (mouse2class2<mTh2_30) & (mouse1class2<mTh1_30)  ) )/length(mouse1class2);
Ath30 = 1/(1-.5)*(fracAgree30-0.5);
mTh1_0 = full(median(mouse1rand));
mTh2_0 = full(median(mouse2rand));
fracAgree0 = ( sum( (mouse2rand>mTh2_0) & (mouse1rand>mTh1_0)  ) ...
    +  sum( (mouse2rand<mTh2_0) & (mouse1rand<mTh1_0)  ) )/length(mouse1rand);
Ath0 = 1/(1-.5)*(fracAgree0-0.5);
plot([mTh1,mTh1],ylim,'k','linewidth',1)
plot(xlim,[mTh2,mTh2],'k','linewidth',1) 
set(gca,'Xtick',full([mTh1,xMax]))
set(gca,'Ytick',full([mTh2,yMax]))
set(gca,'XtickLabel',{'\theta',num2str(round(100*xMax)/100)});%,'fontsize',12)
set(gca,'YtickLabel',{'\theta',num2str(round(100*xMax)/100)});%,'fontsize',12)
ylabel('readout 2','fontsize',12)

display(['untrained A_{0.5}=',num2str(.01*round(100*Ath70)),', ',num2str(.01*round(100*Ath30)),', ',num2str(.01*round(100*Ath0))]);
c1u = corrcoef(mouse2class1,mouse1class1);
c2u = corrcoef(mouse2class2,mouse1class2);
c0u = corrcoef(mouse2rand,mouse1rand);
display(['untrained CC=',num2str(c1u(2)),', ',num2str(c2u(2)),', ',num2str(c0u(2))]);


%% load small values for 4th order response
popDirS = [pwd,'/_simResults/_pirSims/_smallParts/'];
popDLS = dir([popDirS,'pop4thOrderResp*']);
popTotS = load([popDirS,popDLS(1).name]);
for j=2:length(popDLS)
    popTotS(j) = load([popDirS,popDLS(j).name]);
    popTotS(j).rAgr(popTotS(j).rAgr==0)=nan; 
    popTotS(j).rAgr90(popTotS(j).rAgr90==0)=nan; 
    popTotS(j).rAcc(popTotS(j).rAcc==0)=nan; 
    popTotS(j).rCC(popTotS(j).rCC==0)=nan; 
    popTotS(j).rAgrSingle(popTotS(j).rAgrSingle==0)=nan; 
    popTotS(j).rAgrSum(popTotS(j).rAgrSum==0)=nan; 
    popTotS(j).snr(popTotS(j).snr==0)=nan; 
end
smLgCutPt = size(popTotS(1).rAgr,1); %28;

rAgrmeanS = nanmean(popTotS(1).rAgr,3);
for j=2:length(popDLS)
    rAgrmeanS = [rAgrmeanS,nanmean(popTotS(j).rAgr,3)];
end
rAgrmeanS = nanmean(rAgrmeanS,2);

rAgr90meanS = nanmean(popTotS(1).rAgr90,3);
for j=2:length(popDLS)
    rAgr90meanS = [rAgr90meanS,nanmean(popTotS(j).rAgr90,3)];
end
rAgr90meanS = nanmean(rAgr90meanS,2);

rAccmeanS = nanmean(popTotS(1).rAcc,3);
for j=2:length(popDLS)
    rAccmeanS = [rAccmeanS,nanmean(popTotS(j).rAcc,3)];
end
rAccmeanS = nanmean(rAccmeanS,2);

rCCmeanS = nanmean(popTotS(1).rCC,3);
for j=2:length(popDLS)
    rCCmeanS = [rCCmeanS,nanmean(popTotS(j).rCC,3)];
end
rCCmeanS = nanmean(rCCmeanS,2);

rAgrsinglemeanS = nanmean(popTotS(1).rAgrSingle,3);
for j=2:length(popDLS)
    rAgrsinglemeanS = [rAgrsinglemeanS,nanmean(popTotS(j).rAgrSingle,3)];
end
rAgrsinglemeanS = nanmean(rAgrsinglemeanS,2);

rAgrSumS = popTotS(1).rAgrSum;
for j=2:length(popDLS)
    rAgrSumS = [rAgrSumS;popTotS(j).rAgrSum];
end
rAgrSumS = nanmean(rAgrSumS,2);

SNRmeanS = nanmean(popTotS(1).snr,3);
for j=2:length(popDLS)
    SNRmeanS = [SNRmeanS,nanmean(popTotS(j).snr,3)];
end
SNRmeanS = nanmean(SNRmeanS,2);



%% load variables for 4th order response scaling plots

dSteps = popTot(1).dSteps; %[popTotS(1).dSteps, popTot(1).dSteps(2:end)]; %[t1.dSteps,t101.dSteps];

rAgrmean = nanmean(popTot(1).rAgr,3);
for j=2:length(popDL)
    rAgrmean = [rAgrmean,nanmean(popTot(j).rAgr,3)];
end
rAgrmean = [rAgrmeanS(1:smLgCutPt-1); nanmean(rAgrmean(smLgCutPt:end,:),2)];

rAgr90mean = nanmean(popTot(1).rAgr90,3);
for j=2:length(popDL)
    rAgr90mean = [rAgr90mean,nanmean(popTot(j).rAgr90,3)];
end
rAgr90mean = [rAgr90meanS(1:smLgCutPt-1); nanmean(rAgr90mean(smLgCutPt:end,:),2)];
agr90base = popTot(1).agr90base;


rAccmean = nanmean(popTot(1).rAcc,3);
for j=2:length(popDL)
    rAccmean = [rAccmean,nanmean(popTot(j).rAcc,3)];
end
rAccmean = [rAccmeanS(1:smLgCutPt-1); nanmean(rAccmean(smLgCutPt:end,:),2)];


rCCmean = nanmean(popTot(1).rCC,3);
for j=2:length(popDL)
    rCCmean = [rCCmean,nanmean(popTot(j).rCC,3)];
end
rCCmean = [rCCmeanS(1:smLgCutPt-1); nanmean(rCCmean(smLgCutPt:end,:),2)];


SNRmean = nanmean(popTot(1).snr,3);
for j=2:length(popDL)
    SNRmean = [SNRmean,nanmean(popTot(j).snr,3)];
end
SNRmean = [SNRmeanS(1:smLgCutPt-1); nanmean(SNRmean(smLgCutPt:end,:),2)];
SNRraw = SNRmean;
SNRmean = SNRmean - min(SNRmean);
SNRmean = SNRmean/max(SNRmean);
disp(['max SNR = ',num2str(max(SNRraw/99^2))]); % computed SNR with sum rather than mean of 99 test odors, so 99^2 is correction factor

rAgrsinglemean = nanmean(popTot(1).rAgrSingle,3);
for j=2:length(popDL)
    rAgrsinglemean = [rAgrsinglemean,nanmean(popTot(j).rAgrSingle,3)];
end
rAgrsinglemean = [rAgrsinglemeanS(1:smLgCutPt-1); nanmean(rAgrsinglemean(smLgCutPt:end,:),2)];


rAgrSum = popTot(1).rAgrSum;
for j=2:length(popDL)
    rAgrSum = [rAgrSum;popTot(j).rAgrSum];
end
rAgrSum = nanmean(rAgrSum,2);



%try
ZCat = squeeze(popTot(1).ZCat);
for j=2:length(popDL)
    ZCat = [ZCat,squeeze(popTot(j).ZCat)];
end
Z2Cat = squeeze(popTot(1).Z2Cat);
for j=2:length(popDL)
    Z2Cat = [Z2Cat,squeeze(popTot(j).Z2Cat)];
end

%% addional calculations for analytics

accMu = nanmean(squeeze(ZCat),2);
accSigSq = nanvar(squeeze(ZCat),0,2);
SNR = accMu.^2 ./ accSigSq;
%catch
%end

agrFit = zeros(size(rCCmean));
agr90Fit = zeros(size(rCCmean));
readoutAccuracyFit = zeros(size(rCCmean));
%readoutAccuracy = 2*(-.5 + 0.5/innerReps*mean(readoutAccuracy,2)); %NOTE: .5 for z1+z2, 2 for V- and V+

for j=1:length(agrFit)
    try
        agrFit(j) = -1+4*mvncdf([0 0],[0 0],[1 rCCmean(j); rCCmean(j) 1]);
        tmp = mvncdf([1.2815 1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]) + mvncdf([-1.2815 -1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]); % NOTE: normcdf(1.2815) = 0.9
        agr90Fit(j) = (tmp - agr90base)/(1-agr90base); 
    catch
    end
    readoutAccuracyFit(j) = 2*(-0.5 + normcdf(0, -accMu(j), 2*sqrt(accSigSq(j))) );
end
readoutAccuracyFit = nanmean(readoutAccuracyFit, 2);
mA = min(rAccmean);
MA = rAccmean(end)-mA;
mAf = min(readoutAccuracyFit);
MAf = readoutAccuracyFit(end)-mAf;

%% scaling with Ny
lc = [.6,0,.6];
rc = [.6 .3 0];
set(gcf,'defaultAxesColorOrder',[lc; rc]);

axes(bAxes)
semilogx( dSteps, SNRmean, 'linewidth',2, 'color',[.8,.8,.8]);%[.6 .3 0] );
hold on
semilogx( dSteps, rCCmean, '--','linewidth',1, 'color',[.2,.2,.2]);%[.6 .3 0] );
semilogx( dSteps, rAccmean, 'linewidth',1, 'color',[.8 0 .4] );

rAgrmeanNorm = rAgrmean;
rAgrmeanNorm = 2*(rAgrmeanNorm - 0.5);
semilogx( dSteps, rAgrmeanNorm, 'linewidth',1, 'color',[.6 .3 0] );
semilogx( dSteps, 1/(1-agr90base)*(-agr90base+rAgr90mean), '--', 'linewidth',1, 'color',[.9 .6 0] )

agrRef = rAgrmeanNorm;
dStepRef = dSteps;

xlabel('Piriform Inputs per Readout','fontsize',11)
ylabel('Fraction of Maximum','fontsize',11)
[hl,legHndl] = legend('SNR','Correlation','Accuracy','Agreement (\theta=0.5)','Agreement (\theta=0.9)','Location','best'); %'Southeast');
set(hl,'FontSize',10);
legend boxoff
hLines=findobj(legHndl,'type','line');
set(hLines,'linewidth',1.5);
set(gca,'XTick',[1e1,1e2,1e3,1e4,1e5,1e6,1e7])
set(gca,'XTickLabel',{'10^1','','10^3','','10^5','','10^7'})
ylim([0 1])
xlim([400 1e7])
box off


v1 = (dSteps.^2)'./var(ZCat,0,2);
v1 = v1-min(v1);
v1 = v1/max(v1);
figure(paperFig1)


%% agreement scaling loses dependence on theta

sList = linspace(.01,.99999,150);%.1:.1:.9;
thetaList = 1.3;%[.26,.68,1.3];
thetaAdj = normcdf(thetaList);

zeroThreshVal = zeros(length(sList),length(thetaList));
arbitraryThreshVal = zeros(length(sList),length(thetaList));

for i=1:length(sList)
    for j=1:length(thetaList)
        s = sList(i);
        theta = thetaList(j);
        zeroThreshVal(i,j) = -1+4*mvncdf([0 0],[0 0],[1 s; s 1]);
        b = mvncdf([theta theta],[0 0],[1 s; s 1]) + mvncdf([-theta -theta],[0 0],[1 s; s 1]);
        base = thetaAdj(j)^2 + (1-thetaAdj(j))^2;
        arbitraryThreshVal(i,j) = (b-base)/(1-base);
    end
end

muList = .1:.1:2;
ssqList = .1:.1:2;
accGauss = zeros(length(muList),length(ssqList));
snrGauss = zeros(length(muList),length(ssqList));
for i=1:length(muList)
    for j=1:length(ssqList)
        accGauss(i,j) = 2*(-0.5 + normcdf(0, -muList(i), 2*sqrt(ssqList(j))) );
        snrGauss(i,j) = muList(i)^2 / ssqList(j);
    end
end
snrGauss = snrGauss - min(min(snrGauss));
snrGauss = snrGauss/max(max(snrGauss));
accGauss = accGauss - min(min(accGauss));
accGauss = accGauss/max(max(accGauss));

lc = [.6,0,.6];
rc = [.6 .3 0];
set(gcf,'defaultAxesColorOrder',[lc; rc]);

axes(diAxes)
plot(snrGauss, accGauss,'linewidth',1, 'color',[.8 0 .4]);hold on
plot([0 1],[0 1],'k:'); 
box off;
xlabel('SNR', 'fontsize', 11)
ylabel('Accuracy', 'fontsize', 11)
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])


axes(diiAxes)
plot(sList,zeroThreshVal,'linewidth',1, 'color',[.9 .6 0]); hold on
plot(sList,sList,'k:')
box off;
xlabel('Correlation', 'fontsize', 11)
ylabel('Agreement', 'fontsize', 11)
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])


%% FLY scaling with Ny

popDir = [pwd,'/_simResults/_flySims/'];
popDL = dir([popDir,'flyPop4thOrderResp*']);
popTotFly = load([popDir,popDL(1).name]);
for j=2:length(popDL)
    try
        popTotFly(j) = load([popDir,popDL(j).name]);
        popTotFly(j).rAgr(popTotFly(j).rAgr==0)=nan;
        popTotFly(j).rAgr90(popTotFly(j).rAgr90==0)=nan;
        popTotFly(j).rAcc(popTotFly(j).rAcc==0)=nan;
        popTotFly(j).rCC(popTotFly(j).rCC==0)=nan;
        popTotFly(j).rAgrSingle(popTotFly(j).rAgrSingle==0)=nan;
        popTotFly(j).rAgrSum(popTotFly(j).rAgrSum==0)=nan;
        popTotFly(j).snr(popTotFly(j).snr==0)=nan;
    catch
        load([popDir,popDL(j).name]);
        gCC = 0;
        save([popDir,popDL(j).name],'agr90base','T50base','T90base','rAgr90','oCorr','dSteps','kSteps','N4','rAcc','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','REx2','RExU','REx2U','Nsingle','snr');
        popTotFly(j) = load([popDir,popDL(j).name]);
        popTotFly(j).rAgr(popTotFly(j).rAgr==0)=nan;
        popTotFly(j).rAgr90(popTotFly(j).rAgr90==0)=nan;
        popTotFly(j).rAcc(popTotFly(j).rAcc==0)=nan;
        popTotFly(j).rCC(popTotFly(j).rCC==0)=nan;
        popTotFly(j).rAgrSingle(popTotFly(j).rAgrSingle==0)=nan;
        popTotFly(j).rAgrSum(popTotFly(j).rAgrSum==0)=nan;
        popTotFly(j).snr(popTotFly(j).snr==0)=nan;
    end
end

dSteps = popTotFly(1).dSteps; %[popTotS(1).dSteps, popTot(1).dSteps(2:end)]; %[t1.dSteps,t101.dSteps];
rCCmean = nanmean(popTotFly(1).rCC,3);
for j=2:length(popDL)
    rCCmean = [rCCmean,nanmean(popTotFly(j).rCC,3)];
end
rCCmean = nanmean(rCCmean(:,:),2);


snrmean = nanmean(popTotFly(1).snr,3);
for j=2:length(popDL)
    snrmean = [snrmean,nanmean(popTotFly(j).snr,3)];
end
snrmean = nanmean(snrmean(:,:),2);
snrmean = snrmean/max(snrmean);

agrMean50 = nanmean(popTotFly(1).rAgr,3);
for j=2:length(popDL)
    agrMean50 = [agrMean50,nanmean(popTotFly(j).rAgr,3)];
end
agrMean50 = nanmean(agrMean50(:,:),2);

agrMean90 = nanmean(popTotFly(1).rAgr90,3);
for j=2:length(popDL)
    agrMean90 = [agrMean90,nanmean(popTotFly(j).rAgr90,3)];
end
agrMean90 = nanmean(agrMean90(:,:),2);


agrFitFly = zeros(size(agrMean50));
agr90FitFly = zeros(size(agrMean90));

for j=1:length(agrFitFly)
    try
        agrFitFly(j) = -1+4*mvncdf([0 0],[0 0],[1 rCCmean(j); rCCmean(j) 1]);
        tmp = mvncdf([1.2815 1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]) + mvncdf([-1.2815 -1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]); % NOTE: normcdf(1.2815) = 0.9
        agr90FitFly(j) = (tmp - agr90base)/(1-agr90base); 
    catch
    end
end

axes(cAxes)
semilogx( dSteps, rCCmean, 'linewidth',1, 'color',[.2,.2,.2]); hold on
xlabel('Number of Kenyon Cells','fontsize',11)
ylabel('MBON Correlation','fontsize',11)
set(gca,'XTick',[1e1,1e2,1e3,1e4])
set(gca,'XTickLabel',{'10^1','10^2','10^3','10^4'})
set(gca,'YTick',0.3:.1:1)
set(gca,'YTickLabel',{'','0.4','','0.6','','0.8','','1.0'})
ylim([0.3 1])
xl = xlim;
xlim([1e2 2*10^4])
box off

errorbar([333,1666,900],[.559,.842,.758],[.0785,.0402,.0767],'ko','markerfacecolor',[.8 .8 .8],'markersize',8)
text(333*1.2,.59-.04,'\alpha''2','color',[0 0 0])
text(1666*1.2,.88-.04,'\gamma1pedc','color',[0 0 0])
text(1000*1.2,.8-.04,'\alpha2sc','color',[0 0 0])

[hl,legHndl] = legend('Model','Hige et al.','Location','Southeast');
set(hl,'FontSize',10);
legend boxoff
hLines=findobj(legHndl,'type','line');
set(hLines,'linewidth',1.5);


