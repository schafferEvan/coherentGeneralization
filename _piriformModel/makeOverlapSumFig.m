

function makeOverlapSumFig
global paperFig classLevel1 classLevel2 class1Color class2Color

% updated to use only "correlated" odors instead of mix & class on 12/29/15
% split into plotting part and data generating part (calculateSumFigParts.m) on 1/5/16

%load('partsForSumFig_Rand.mat');
%load('partsForSumFig.mat');

%load('partsForSumFig_100kNy.mat');
load('overlapForSumFig.mat');

glomOverlap = glomOverlap(~isnan(glomOverlap));
pirOverlap = pirOverlap(~isnan(pirOverlap));


%% FIGURE 5 ---------------------------------------------------------



paperFig = figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=600;%pos(3)*1.5;
pos(4)=600;%pos(4)*1.5;
set(gcf,'Position',pos)
%blueMap = [0*(1:64)',1/(3*64)*(64:-1:1)',(1/64)*(64:-1:1)'];
%greenMap = [0*(1:64)',(1/64)*(64:-1:1)',(1/64)*(1:64)'];
%colormap([greenMap;blueMap;hot(64)])




axes('position',[.15  .43  .4  .1]); hold on
hrange = 0:max(glomOverlap)/100:max(glomOverlap); %0:.01:1; %hrange = 0:1:201; 

hout = hist(glomOverlap(corrIsRandIdx),hrange)/sum(corrIsRandIdx);
hin1 = hist(glomOverlap(corrIsClass1Idx),hrange)/sum(corrIsClass1Idx);
hin2 = hist(glomOverlap(corrIsClass2Idx),hrange)/sum(corrIsClass2Idx);
%hmix = hist(glomOverlap(corrIsMixtureIdx),hrange)/sum(corrIsMixtureIdx);
% hout = hist(glomRepSimilarity(similarityIsOtherIdx),hrange)/sum(similarityIsOtherIdx);
% hin = hist(glomRepSimilarity(similarityIsClassIdx),hrange)/sum(similarityIsClassIdx);
% hmix = hist(glomRepSimilarity(similarityIsMixtureIdx),hrange)/sum(similarityIsMixtureIdx);
% %fill([hrange(1),hrange,hrange(end)],-hx,'color',[0,1,1],'linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hout,0],'k','linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hin2,0],class2Color,'linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hin1,0],class1Color,'linewidth',1); hold on
%fill([hrange(1),hrange,hrange(end)],[0,-hin1,0],[0 classLevel1+greenColorAdjust 0],'linewidth',1); hold on
%fill([hrange(1),hrange,hrange(end)],[0,-hmix,0], [0 .333 1],'linewidth',1); hold on
%alpha(.7); % (transparency)
xlim([0 0.08]);
%xlim([0, max(hrange)])
ylim([-.2 0])
box off;
axis off
%xlabel('Overlap in Representation','fontsize',12)
ylabel('P(r)','fontsize',12)
%title('Odor Similarity (raw correlation)', 'fontsize',12);



axes('position',[.05  .53  .1  .4]); hold on
%subplot(3,6,6+(5:6));
hrange = 0:max(pirOverlap)/100:max(pirOverlap); %0:.01:1; %hrange = 0:1:201; 

hout = hist(pirOverlap(corrIsRandIdx),hrange)/sum(corrIsRandIdx);
hin1 = hist(pirOverlap(corrIsClass1Idx),hrange)/sum(corrIsClass1Idx);
hin2 = hist(pirOverlap(corrIsClass2Idx),hrange)/sum(corrIsClass2Idx);
%hmix = hist(pirOverlap(corrIsMixtureIdx),hrange)/sum(corrIsMixtureIdx);
% hx = hist(pirRepSimilarity(similarityIsCross),hrange)/sum(similarityIsCross);
% hout = hist(pirRepSimilarity(similarityIsOtherIdx),hrange)/sum(similarityIsOtherIdx);
% hin = hist(pirRepSimilarity(similarityIsClassIdx),hrange)/sum(similarityIsClassIdx);
% hmix = hist(pirRepSimilarity(similarityIsMixtureIdx),hrange)/sum(similarityIsMixtureIdx);
% %fill(-hx,[0,hrange,0],'color',[0,1,1],'linewidth',1); hold on
fill([0,0,-hout,0,0],[0,hrange(1),hrange,hrange(end),0],'k','linewidth',1); hold on
fill([0,-hin2,0],[hrange(1),hrange,hrange(end)],class2Color,'linewidth',1); hold on
fill([0,-hin1,0],[hrange(1),hrange,hrange(end)],class1Color,'linewidth',1); hold on
%fill([0,-hmix,0],[hrange(1),hrange,hrange(end)],[0 .333 1],'linewidth',1); hold on
%alpha(.7); % (transparency)
ylim([0 0.08]);
%ylim([0, max(hrange)])
xlim([-.2 0])
box off;
axis off
%ylabel('Overlap in Representation','fontsize',12)
xlabel('P(r)','fontsize',12)
%title('Odor Similarity (raw correlation)', 'fontsize',12);



% correlations -------------
axes('position',[.15  .53  .4  .4]); hold on

%corrIsMixtureLoc = find(corrIsMixtureIdx);
corrIsClass1Loc = find(corrIsClass1Idx);
corrIsClass2Loc = find(corrIsClass2Idx);
corrIsOtherLoc = find(corrIsRandIdx);

% for j=1:npieces
%     thisPiece = ( mixCorrSimilarity > (j-1)/npieces ) & ( mixCorrSimilarity < j/npieces );
%     %plot( glomOverlap( thisPiece ) ,    pirOverlap( thisPiece ),'.', 'color', [.75*j/npieces .5*j/npieces j/npieces],'markersize',5); hold on
%     plot( glomOverlap( thisPiece ) ,    pirOverlap( thisPiece ),'.', 'color', [0 0 .5+.5*tanh(-3+6*j/npieces)],'markersize',5); hold on
% end

plot( glomOverlap( corrIsMixtureIdx ) ,    pirOverlap( corrIsMixtureIdx ),'.', 'color', [.87 .87 .87],'markersize',12); hold on
plot(glomOverlap(  corrIsClass1Loc ) ,  pirOverlap( corrIsClass1Loc ),'.', 'color', class1Color,'markersize',5)
plot(glomOverlap(  corrIsClass2Loc ) ,  pirOverlap( corrIsClass2Loc ),'.', 'color', class2Color,'markersize',5)
plot(glomOverlap(  corrIsOtherLoc ) ,  pirOverlap( corrIsOtherLoc ),'k.','markersize',5)

box off
xlabel({'';'Overlap in glomerulus representation'},'fontsize',12)
ylabel({'Overlap in piriform representation';''},'fontsize',12)
xlim([0 0.08]);
ylim([0 0.08]);



%% FIGURE 5 ---------------------------------------------------------



paperFig = figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=600;%pos(3)*1.5;
pos(4)=600;%pos(4)*1.5;
set(gcf,'Position',pos)
%blueMap = [0*(1:64)',1/(3*64)*(64:-1:1)',(1/64)*(64:-1:1)'];
%greenMap = [0*(1:64)',(1/64)*(64:-1:1)',(1/64)*(1:64)'];
%colormap([greenMap;blueMap;hot(64)])




axes('position',[.15  .43  .4  .1]); hold on
hrange = 0:max(glomCorr)/100:max(glomCorr); %0:.01:1; %hrange = 0:1:201; 

hout = hist(glomCorr(corrIsRandIdx),hrange)/sum(corrIsRandIdx);
hin1 = hist(glomCorr(corrIsClass1Idx),hrange)/sum(corrIsClass1Idx);
hin2 = hist(glomCorr(corrIsClass2Idx),hrange)/sum(corrIsClass2Idx);
%hmix = hist(glomCorr(corrIsMixtureIdx),hrange)/sum(corrIsMixtureIdx);
% hout = hist(glomRepSimilarity(similarityIsOtherIdx),hrange)/sum(similarityIsOtherIdx);
% hin = hist(glomRepSimilarity(similarityIsClassIdx),hrange)/sum(similarityIsClassIdx);
% hmix = hist(glomRepSimilarity(similarityIsMixtureIdx),hrange)/sum(similarityIsMixtureIdx);
% %fill([hrange(1),hrange,hrange(end)],-hx,'color',[0,1,1],'linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hout,0],'k','linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hin2,0],class2Color,'linewidth',1); hold on
fill([hrange(1),hrange,hrange(end)],[0,-hin1,0],class1Color,'linewidth',1); hold on
%fill([hrange(1),hrange,hrange(end)],[0,-hin1,0],[0 classLevel1+greenColorAdjust 0],'linewidth',1); hold on
%fill([hrange(1),hrange,hrange(end)],[0,-hmix,0], [0 .333 1],'linewidth',1); hold on
%alpha(.7); % (transparency)
xlim([0 0.08]);
%xlim([0, max(hrange)])
ylim([-.2 0])
box off;
axis off
%xlabel('Overlap in Representation','fontsize',12)
ylabel('P(r)','fontsize',12)
%title('Odor Similarity (raw correlation)', 'fontsize',12);



axes('position',[.05  .53  .1  .4]); hold on
%subplot(3,6,6+(5:6));
hrange = 0:max(pirOverlap)/100:max(pirOverlap); %0:.01:1; %hrange = 0:1:201; 

hout = hist(pirOverlap(corrIsRandIdx),hrange)/sum(corrIsRandIdx);
hin1 = hist(pirOverlap(corrIsClass1Idx),hrange)/sum(corrIsClass1Idx);
hin2 = hist(pirOverlap(corrIsClass2Idx),hrange)/sum(corrIsClass2Idx);
%hmix = hist(pirOverlap(corrIsMixtureIdx),hrange)/sum(corrIsMixtureIdx);
% hx = hist(pirRepSimilarity(similarityIsCross),hrange)/sum(similarityIsCross);
% hout = hist(pirRepSimilarity(similarityIsOtherIdx),hrange)/sum(similarityIsOtherIdx);
% hin = hist(pirRepSimilarity(similarityIsClassIdx),hrange)/sum(similarityIsClassIdx);
% hmix = hist(pirRepSimilarity(similarityIsMixtureIdx),hrange)/sum(similarityIsMixtureIdx);
% %fill(-hx,[0,hrange,0],'color',[0,1,1],'linewidth',1); hold on
fill([0,0,-hout,0,0],[0,hrange(1),hrange,hrange(end),0],'k','linewidth',1); hold on
fill([0,-hin2,0],[hrange(1),hrange,hrange(end)],class2Color,'linewidth',1); hold on
fill([0,-hin1,0],[hrange(1),hrange,hrange(end)],class1Color,'linewidth',1); hold on
%fill([0,-hmix,0],[hrange(1),hrange,hrange(end)],[0 .333 1],'linewidth',1); hold on
%alpha(.7); % (transparency)
ylim([0 0.08]);
%ylim([0, max(hrange)])
xlim([-.2 0])
box off;
axis off
%ylabel('Overlap in Representation','fontsize',12)
xlabel('P(r)','fontsize',12)
%title('Odor Similarity (raw correlation)', 'fontsize',12);



% correlations -------------
axes('position',[.15  .53  .4  .4]); hold on

%corrIsMixtureLoc = find(corrIsMixtureIdx);
corrIsClass1Loc = find(corrIsClass1Idx);
corrIsClass2Loc = find(corrIsClass2Idx);
corrIsOtherLoc = find(corrIsRandIdx);

% for j=1:npieces
%     thisPiece = ( mixCorrSimilarity > (j-1)/npieces ) & ( mixCorrSimilarity < j/npieces );
%     %plot( glomCorr( thisPiece ) ,    pirOverlap( thisPiece ),'.', 'color', [.75*j/npieces .5*j/npieces j/npieces],'markersize',5); hold on
%     plot( glomCorr( thisPiece ) ,    pirOverlap( thisPiece ),'.', 'color', [0 0 .5+.5*tanh(-3+6*j/npieces)],'markersize',5); hold on
% end

plot( glomCorr( corrIsMixtureIdx ) ,    pirOverlap( corrIsMixtureIdx ),'.', 'color', [.87 .87 .87],'markersize',12); hold on
plot(glomCorr(  corrIsClass1Loc ) ,  pirOverlap( corrIsClass1Loc ),'.', 'color', class1Color,'markersize',5)
plot(glomCorr(  corrIsClass2Loc ) ,  pirOverlap( corrIsClass2Loc ),'.', 'color', class2Color,'markersize',5)
plot(glomCorr(  corrIsOtherLoc ) ,  pirOverlap( corrIsOtherLoc ),'k.','markersize',5)

box off
xlabel({'';'Correlation in glomerulus representation'},'fontsize',12)
ylabel({'Overlap in piriform representation';''},'fontsize',12)
xlim([0 1]);
ylim([0 0.08]);

