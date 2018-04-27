

function makeFlyScalingFig
global paperFig class1Color
initializeGlobals_Server;


%% load variables for 4th order response scaling plots
popDir = [pwd,'/_newBig0_aHebb/'];%.7/'];%'~/Dropbox/_AxelLab/matlab/_dan/_pop4thOrderResp_0.7/';
%popDir = [pwd,'/_Sc=9/'];%'~/Dropbox/_AxelLab/matlab/_dan/_pop4thOrderResp_0.7/';
popDL = dir([popDir,'flyPop4thOrderResp*']);
popTot = load([popDir,popDL(1).name]);
for j=2:length(popDL)
    try
        popTot(j) = load([popDir,popDL(j).name]);
        popTot(j).rAgr(popTot(j).rAgr==0)=nan;
        popTot(j).rAgr90(popTot(j).rAgr90==0)=nan;
        popTot(j).rAcc(popTot(j).rAcc==0)=nan;
        popTot(j).rCC(popTot(j).rCC==0)=nan;
        popTot(j).rAgrSingle(popTot(j).rAgrSingle==0)=nan;
        popTot(j).rAgrSum(popTot(j).rAgrSum==0)=nan;
        popTot(j).snr(popTot(j).snr==0)=nan;
        %popTot(j).ZCat(popTot(j).ZCat==0)=nan;
        %popTot(j).Z2Cat(popTot(j).Z2Cat==0)=nan;
    catch
        load([popDir,popDL(j).name]);
        gCC = 0;
        save([popDir,popDL(j).name],'agr90base','T50base','T90base','rAgr90','oCorr','dSteps','kSteps','N4','rAcc','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','REx2','RExU','REx2U','Nsingle','snr');
        popTot(j) = load([popDir,popDL(j).name]);
        popTot(j).rAgr(popTot(j).rAgr==0)=nan;
        popTot(j).rAgr90(popTot(j).rAgr90==0)=nan;
        popTot(j).rAcc(popTot(j).rAcc==0)=nan;
        popTot(j).rCC(popTot(j).rCC==0)=nan;
        popTot(j).rAgrSingle(popTot(j).rAgrSingle==0)=nan;
        popTot(j).rAgrSum(popTot(j).rAgrSum==0)=nan;
        popTot(j).snr(popTot(j).snr==0)=nan;
        %popTot(j).ZCat(popTot(j).ZCat==0)=nan;
        %popTot(j).Z2Cat(popTot(j).Z2Cat==0)=nan;
    end
end


dSteps = popTot(1).dSteps; %[popTotS(1).dSteps, popTot(1).dSteps(2:end)]; %[t1.dSteps,t101.dSteps];

rAgrmean = nanmean(popTot(1).rAgr,3);
for j=2:length(popDL)
    rAgrmean = [rAgrmean,nanmean(popTot(j).rAgr,3)];
end
%rAgrmean = (99/98)*nanmean(rAgrmean(:,:),2);
rAgrmean = nanmean(rAgrmean(:,:),2);
% NOTE: correction factor of 99/98 is because median point is always on
% boundary and thus technically never in "agreement".

rAgr90mean = nanmean(popTot(1).rAgr90,3);
for j=2:length(popDL)
    rAgr90mean = [rAgr90mean,nanmean(popTot(j).rAgr90,3)];
end
rAgr90mean = nanmean(rAgr90mean(:,:),2);
agr90base = popTot(1).agr90base;


rAccmean = nanmean(popTot(1).rAcc,3);
for j=2:length(popDL)
    rAccmean = [rAccmean,nanmean(popTot(j).rAcc,3)];
end
rAccmean = nanmean(rAccmean(:,:),2);


rCCmean = nanmean(popTot(1).rCC,3);
for j=2:length(popDL)
    rCCmean = [rCCmean,nanmean(popTot(j).rCC,3)];
end
rCCmean = nanmean(rCCmean(:,:),2);


SNRmean = nanmean(popTot(1).snr,3);
for j=2:length(popDL)
    SNRmean = [SNRmean,nanmean(popTot(j).snr,3)];
end
SNRmean = nanmean(SNRmean(:,:),2);
%SNRmean = SNRmean - min(SNRmean);
SNRmean = SNRmean/max(SNRmean);


rAgrsinglemean = nanmean(popTot(1).rAgrSingle,3);
for j=2:length(popDL)
    rAgrsinglemean = [rAgrsinglemean,nanmean(popTot(j).rAgrSingle,3)];
end
rAgrsinglemean = nanmean(rAgrsinglemean(:,:),2);


rAgrSum = popTot(1).rAgrSum;
for j=2:length(popDL)
    rAgrSum = [rAgrSum;popTot(j).rAgrSum];
end
rAgrSum = nanmean(rAgrSum,2);



%try
% ZCat = squeeze(popTot(1).ZCat);
% for j=2:length(popDL)
%     ZCat = [ZCat,squeeze(popTot(j).ZCat)];
% end
% Z2Cat = squeeze(popTot(1).Z2Cat);
% for j=2:length(popDL)
%     Z2Cat = [Z2Cat,squeeze(popTot(j).Z2Cat)];
% end

%ZCat = [ZCatS; [ZCat(2:end,:), nan(size(ZCat,1)-1, size(ZCatS,2)-size(ZCat,2))  ] ];
%Z2Cat = [Z2CatS; [Z2Cat(2:end,:), nan(size(Z2Cat,1)-1, size(Z2CatS,2)-size(Z2Cat,2))  ] ];

%% addional calculations for analytics

%accMu = nanmean(squeeze(ZCat),2);
%accSigSq = nanvar(squeeze(ZCat),0,2);
%SNR = accMu.^2 ./ accSigSq;
%catch
%end

% agrFit = zeros(size(rCCmean));
% agr90Fit = zeros(size(rCCmean));
% readoutAccuracyFit = zeros(size(rCCmean));
% %readoutAccuracy = 2*(-.5 + 0.5/innerReps*mean(readoutAccuracy,2)); %NOTE: .5 for z1+z2, 2 for V- and V+
% 
% for j=1:length(agrFit)
%     try
%         agrFit(j) = -1+4*mvncdf([0 0],[0 0],[1 rCCmean(j); rCCmean(j) 1]);
%         tmp = mvncdf([1.2815 1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]) + mvncdf([-1.2815 -1.2815],[0 0],[1 rCCmean(j); rCCmean(j) 1]); % NOTE: normcdf(1.2815) = 0.9
%         agr90Fit(j) = (tmp - agr90base)/(1-agr90base); 
%     catch
%     end
%     %readoutAccuracyFit(j) = 2*(-0.5 + normcdf(0, -accMu(j), 2*sqrt(accSigSq(j))) );
% end
% readoutAccuracyFit = nanmean(readoutAccuracyFit, 2);
% mA = min(rAccmean);
% MA = rAccmean(end)-mA;
% mAf = min(readoutAccuracyFit);
% MAf = readoutAccuracyFit(end)-mAf;




%% 1st Model Figure


figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=800;%pos(3)*1.5;
pos(4)=400;%pos(4)*1.5;
set(gcf,'Position',pos)



for j=1:2
    axes('position', [.05+(j-1)*.5  .15  .4  .8] ); hold on
    
    if j==1
        mouse1class1 = popTot(5).RExU(end,:);
        mouse2class1 = popTot(5).REx2U(end,:);
    else
        mouse1class1 = popTot(5).REx(end,:);
        mouse2class1 = popTot(5).REx2(end,:);
    end
    
    xMax = full(max(mouse1class1));
    yMax = full(max(mouse2class1));
    xMin = full(min(mouse1class1));
    yMin = full(min(mouse2class1));
    
    mTh1 = full(median(mouse1class1));
    mTh2 = full(median(mouse2class1));
    rectangle('Position',[mTh1,yMin,xMax-mTh1,mTh2-yMin],'facecolor',.9*[1,1,1],'edgecolor','none')
    rectangle('Position',[xMin,mTh2,mTh1-xMin,yMax-mTh2],'facecolor',.9*[1,1,1],'edgecolor','none')
    
    plot(mouse1class1,mouse2class1,'.','color',class1Color,'markersize',6); %[.7 .7 .7]
    box off
    ylim([yMin yMax]);
    xlim([xMin xMax]);
    fracAgree = ( sum( (mouse2class1>mTh2) & (mouse1class1>mTh1)  ) ...
        +  sum( (mouse2class1<mTh2) & (mouse1class1<mTh1)  ) )/length(mouse1class1);
    plot([mTh1,mTh1],ylim,'k','linewidth',1)
    plot(xlim,[mTh2,mTh2],'k','linewidth',1)
    set(gca,'Xtick',full([mTh1,xMax]))
    set(gca,'Ytick',full([mTh2,yMax]))
    set(gca,'XtickLabel',{'',''},'fontsize',12)
    set(gca,'YtickLabel',{'\theta','1'},'fontsize',12)
    ylabel('MBON 2','fontsize',12)
    xlabel('MBON 1','fontsize',12)
    set(gca,'XtickLabel',{'\theta','1'},'fontsize',12)
    %text( .55*(xMax-xMin)+xMin, .21*(yMax-yMin)+yMin, ['N_y=',NexListStr{j}],'color',[0 0 0]);
    text( .55*(xMax-xMin)+xMin, .07*(yMax-yMin)+yMin, ['A=',num2str(.01*round(100*fracAgree))],'color',[.6 .3 0]);
    %text( .03*(xMax-xMin)+xMin, .91*(yMax-yMin)+yMin, ['N_y=',NexListStr{j}],'color',[0 0 0]);
    %text( .03*(xMax-xMin)+xMin, .75*(yMax-yMin)+yMin, ['A=',num2str(.01*round(100*fracAgree))],'color',[.6 .3 0]);
    if j==1
        title('Untrained')
    else
        title('Trained')
    end
end



%% 2nd Model Figure


paperFig = figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(1)=1000;
pos(3)=900;%pos(3)*1.5;
pos(4)=600;%pos(4)*1.5;
set(gcf,'Position',pos)


NexList = dSteps([24,28,30,33]);%[1e3,1e4,1e5];
NexListStr = {num2str(NexList(1)),num2str(NexList(2)),num2str(NexList(3)),num2str(NexList(4))}; %{'10^2','10^3','10^4'};
exIdx = [find(popTot(1).dSteps==NexList(1)),find(popTot(1).dSteps==NexList(2)),find(popTot(1).dSteps==NexList(3)),find(popTot(1).dSteps==NexList(4))];
xSz = 0.15;
axesList = [.15  .78  xSz*pos(4)/pos(3)  xSz;...
            .15  .61  xSz*pos(4)/pos(3)  xSz;...
            .15  .44  xSz*pos(4)/pos(3)  xSz;...
            .15  .27  xSz*pos(4)/pos(3)  xSz];
        
for j=1:length(exIdx)
    axes('position',axesList(j,:)); hold on
    
    mouse1class1 = popTot(5).REx(exIdx(j),:);
    mouse2class1 = popTot(5).REx2(exIdx(j),:);
    
    xMax = full(max(mouse1class1));
    yMax = full(max(mouse2class1));
    xMin = full(min(mouse1class1));
    yMin = full(min(mouse2class1));
    
    mTh1 = full(median(mouse1class1));
    mTh2 = full(median(mouse2class1));
    rectangle('Position',[mTh1,yMin,xMax-mTh1,mTh2-yMin],'facecolor',.9*[1,1,1],'edgecolor','none')
    rectangle('Position',[xMin,mTh2,mTh1-xMin,yMax-mTh2],'facecolor',.9*[1,1,1],'edgecolor','none')
    
    plot(mouse1class1,mouse2class1,'.','color',class1Color,'markersize',6); %[.7 .7 .7]
    box off
    ylim([yMin yMax]);
    xlim([xMin xMax]);
    fracAgree = ( sum( (mouse2class1>mTh2) & (mouse1class1>mTh1)  ) ...
        +  sum( (mouse2class1<mTh2) & (mouse1class1<mTh1)  ) )/length(mouse1class1);
    plot([mTh1,mTh1],ylim,'k','linewidth',1)
    plot(xlim,[mTh2,mTh2],'k','linewidth',1)
    set(gca,'Xtick',full([mTh1,xMax]))
    set(gca,'Ytick',full([mTh2,yMax]))
    set(gca,'XtickLabel',{'',''},'fontsize',12)
    set(gca,'YtickLabel',{'\theta','1'},'fontsize',12)
    ylabel('MBON 2','fontsize',12)
    text( .55*(xMax-xMin)+xMin, .21*(yMax-yMin)+yMin, ['N_y=',NexListStr{j}],'color',[0 0 0]);
    text( .55*(xMax-xMin)+xMin, .07*(yMax-yMin)+yMin, ['A=',num2str(.01*round(100*fracAgree))],'color',[.6 .3 0]);
    %text( .03*(xMax-xMin)+xMin, .91*(yMax-yMin)+yMin, ['N_y=',NexListStr{j}],'color',[0 0 0]);
    %text( .03*(xMax-xMin)+xMin, .75*(yMax-yMin)+yMin, ['A=',num2str(.01*round(100*fracAgree))],'color',[.6 .3 0]);
end
xlabel('MBON 1','fontsize',12)
set(gca,'XtickLabel',{'\theta','1'},'fontsize',12)





%% scaling with Ny
lc = [.6,0,.6];
rc = [.6 .3 0];
set(gcf,'defaultAxesColorOrder',[lc; rc]);

axes('position',[.36  .27  .52  .66]);

semilogx( dSteps, SNRmean, 'linewidth',2, 'color',[.8,.8,.8]); hold on
%try
%SNR = SNR - min(SNR);
%SNR = SNR/max(SNR);
%semilogx( dSteps*10, SNR, ':','linewidth',1, 'color',[0,0,0]);%[.6 .3 0] );
%semilogx( dSteps*10, readoutAccuracyFit, '--', 'linewidth',1, 'color',[.8 0 .4] ); % NOTE: this may work better with calc based on individual reps, instead of post hoc, as is the case with snr
%catch
%end
semilogx( dSteps, rCCmean, 'linewidth',1, 'color',[.2,.2,.2]); hold on
%semilogx( dSteps, rAccmean, 'linewidth',1, 'color',[.8 0 .4] );

rAgrmeanNorm = rAgrmean;
rAgrmeanNorm = 2*(rAgrmeanNorm - 0.5);
semilogx( dSteps, rAgrmeanNorm, 'linewidth',1, 'color',[.6 .3 0] );
semilogx( dSteps, 1/(1-agr90base)*(-agr90base+rAgr90mean), 'linewidth',1, 'color',[.9 .6 0] )
%semilogx( dSteps, agrFit, '--', 'linewidth',1, 'color',[.6 .3 0] )
%semilogx( dSteps, agr90Fit, '--', 'linewidth',1, 'color',[.9 .6 0] )

xlabel('Number of Kenyon Cells','fontsize',12)
ylabel('Correlation between two MBONs','fontsize',12)
%[hl,legHndl] = legend('SNR','Correlation','Accuracy','Agreement (\theta=0.5)','Agreement (\theta=0.9)','Location','Southeast');
set(gca,'XTick',[1e1,1e2,1e3,1e4])
set(gca,'XTickLabel',{'10^1','10^2','10^3','10^4'})
set(gca,'YTick',0:.1:1)
set(gca,'YTickLabel',{'','','0.2','','0.4','','0.6','','0.8','','1.0'})
ylim([0 1])
xl = xlim;
xlim([1e1 2*10^4])
box off

%errorbar([333,1666,1000],[.59,.88,.85],[.15,.08,.15],'ko','markerfacecolor',[0 0 1],'markersize',10)
errorbar([333,1666,1000],[.59,.88,.8],[.15,.08,.1],'ko','markerfacecolor',[0 0 1],'markersize',10)
text(333*1.08,.59-.03,'\alpha''2','color','b')
text(1666*1.08,.88-.03,'\gamma1pedc','color','b')
text(1000*1.08,.8-.03,'\alpha2sc','color','b')
set(gca,'fontsize',14)

%[hl,legHndl] = legend('SNR','Correlation','Agreement','','Hige et al., Nature, 2015','Location','Southeast');
[hl,legHndl] = legend('Model Correlation','Hige et al., Nature, 2015','Location','Southeast');
set(hl,'FontSize',14);
legend boxoff
hLines=findobj(legHndl,'type','line');
set(hLines,'linewidth',2);


