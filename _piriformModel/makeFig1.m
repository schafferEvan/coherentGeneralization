
function makeFig1
global Ny Sp Sx figureUno classLevel1 classLevel2 greenColorAdjust class1Color class2Color

% switched to corr odors and removed mix on 12/25/15
% added variability to input sparseness for comparison to data on 12/24/17

initializeGlobals;
tmpTime = clock;
rng(1e5*tmpTime(end)); %seed random number generator


No = 100;
figureUno = figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
%pos(3)=600;%pos(3)*1.5;
pos(4) = pos(4)*2; %pos(4)=600;%
set(gcf,'Position',pos)
blueMap = [0*(1:64)',1/(1.5*64)*(64:-1:1)',(1/64)*(64:-1:1)'];
%greenMap = [0*(1:64)',(1/64)*(1:64)',0*(1:64)'];
greenMap = [repmat(class2Color,32,1);repmat(class1Color,32,1)];
colormap([greenMap;blueMap;hot(64)])
%greenMap = [0*(1:64)',.8*(1/64)*(64:-1:1)',0*(1:64)'];



%% data panel
%[nCellsA,nCellsB,nCellsTot,observedOverlap,expectedOverlap,compNames] = plotDataOverlap;
%display([num2str(observedOverlap(~isnan(observedOverlap))),repmat('  ',length(observedOverlap(~isnan(observedOverlap))),1),num2str(expectedOverlap(~isnan(expectedOverlap))),repmat('  ',length(observedOverlap(~isnan(observedOverlap))),1),char(compNames(~isnan(expectedOverlap)))]);

xlFile = '/Users/evan/Dropbox/_AxelLab/matlab/_dan/_data/pairwise_comparison_b.xls';
sData = xlsread(xlFile);
%figure; plot(sData(34:47,end-1),sData(34:47,end),'b.');
%hold all; plot(sData(1:24,end-1),sData(1:24,end),'r.')
expectedOverlapMix = sData(34:47,end-1);
observedOverlapMix = sData(34:47,end);
expectedOverlapMono = sData(1:24,end-1);
observedOverlapMono = sData(1:24,end);
monoLabels = {'oct vs cit'; 'oct vs ems'; 'hex vs oct'};
mixLabels = {'oct vs ems+oct'; 'oct vs eug+oct'; 'oct vs ace+oct'};
mixLabelIdx = [5,10,1];


figure(figureUno);
annotation('textbox', [.1  .55  .35  .02],...
    'String', 'Data','BackgroundColor','none','Color','k',...
    'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');
annotation('textbox', [.55  .55  .35  .02],...
    'String', 'Model','BackgroundColor','none','Color','k',...
    'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');

%subplot(3,12,[13:16,25:28]); hold on
axes('position',[.1  .33  .35  .2]); hold on
plot(expectedOverlapMono,observedOverlapMono,'.','markersize',20,'color',[0 0 0]);
for j=1:3
    text(expectedOverlapMono(j)-.00035-.003*(j==3), observedOverlapMono(j)-.001*(-1)^j, monoLabels{j});
end
box off;
%xlabel('Chance Frequency','fontsize',12)
ylabel('Observed Frequency','fontsize',12)
xlim([0 .02]);
ylim([0 .03]);
set(gca,'XTick',0:.01:.02,'fontsize',12);
set(gca,'YTick',0:.01:.03,'fontsize',12);
%set(gca,'YTick',0:.02:.1,'fontsize',12);
lineTop = min( max(xlim),max(ylim) );
plot([0 lineTop],[0 lineTop],'color',[.5 .5 .5],'linewidth',1)


axes('position',[.1  .05  .35  .2]); hold on
plot(expectedOverlapMix,observedOverlapMix,'.','markersize',20,'color',[0 0 0]);
for j=1:3
    text(expectedOverlapMix(mixLabelIdx(j))-.00035-.007*(j==3), observedOverlapMix(mixLabelIdx(j))-.0065*(-1)^j, mixLabels{j});
end
box off;
xlabel('Chance Frequency','fontsize',12)
ylabel('Observed Frequency','fontsize',12)
xlim([0 .02]);
ylim([0 .12]);
set(gca,'XTick',0:.01:.02,'fontsize',12);
set(gca,'YTick',0:.05:.15,'fontsize',12);
lineTop = min( max(xlim),max(ylim) );
plot([0 lineTop],[0 lineTop],'color',[.5 .5 .5],'linewidth',1)


%% model


% first make correlation matrix:
NoTot = 3*No;
Sclass = classLevel2*(ones(2*No)-eye(2*No))+eye(2*No);
Sclass(No+1:end,No+1:end) = classLevel1*(ones(No)-eye(No))+eye(No);
S = eye(NoTot);
S(No+1:end,No+1:end) = Sclass;

[glomRepTemp,odorSimilarity] = makeOdors('correlated',S,NoTot);
glomerulusRepRand = glomRepTemp(:,1:No);
similarityRand = -eps*(odorSimilarity(1:No)<0); 
glomerulusRepCorr2 = glomRepTemp(:,No+1:2*No);
similarityClass2 = -2 + odorSimilarity(No+1:2*No);
glomerulusRepCorr1 = glomRepTemp(:,2*No+1:end);
similarityClass1 = -2 + odorSimilarity(2*No+1:end);

% %random
% [glomerulusRepRand,similarityRand] = makeOdors('classWithoutTemplate',0,No);
% similarityRand = -eps*(similarityRand<0); %-eps*(similarityRand<0);
[piriformRepRand,~,~,J] = makePiriform(glomerulusRepRand);

%graded
%[glomerulusRepCorr1,similarityClass1] = makeOdors('correlated',classLevel1,No); %classWithoutTemplate
%similarityClass1 = -2 + similarityClass1 + greenColorAdjust;
piriformRepClass1 = makePiriform(glomerulusRepCorr1,Sp,J);
%similarityClass1 = -1+similarityClass1; %-2+eps*(similarityClass1<0);

%[glomerulusRepCorr2,similarityClass2] = makeOdors('correlated',classLevel2,No); %classWithoutTemplate
%similarityClass2 = -2 + similarityClass2 + greenColorAdjust;
piriformRepClass2 = makePiriform(glomerulusRepCorr2,Sp,J);
%similarityClass2 = -1+similarityClass2; %-2+eps*(similarityClass1<0);

%[glomRepTempMix,odorSimilarity] = makeOdors('mix',0,No);
%mixtureOdorSimilarity = -1+odorSimilarity/No;
%piriformRepMix = makePiriform(glomRepTempMix,Sp,J);



%combined for figure
glomerulusRep = [glomerulusRepRand,glomerulusRepCorr2,glomerulusRepCorr1];%,glomRepTempMix];
piriformRep = [piriformRepRand,piriformRepClass2,piriformRepClass1];%,piriformRepMix];

z = load('partsForFig1odorPairComparison.mat');

%% FIGURE 1G model ---------------------------------------



% FIGURE 1 REPRESENTATIONS -----------------------------------------

figure(figureUno)
oPlotIdx = 1:3*No; %[1:50,101:150,201:250,301:350];

% glomeruli ---------
gPlotIdx = [1:15, 31:50, 101:115];%1:150; %75:225;
glomerulusRepScaledForFig = glomerulusRep(gPlotIdx,oPlotIdx)/max(.7*max(glomerulusRep(gPlotIdx,oPlotIdx)));
%axes('position',[.4  .7  .25  .25])
axes('position',[.25  .78  .5  .15])
%subplot(3,12,(8:9));
%subplot(3,6,(2:3));
imagesc(glomerulusRepScaledForFig);
box off;
%xlabel('Odor','fontsize',12)
ylabel('Glomerulus','fontsize',12); % # (1000 Total)
%title('Bulb','fontsize',12); %Glomerular Representation
set(gca,'YTick',[])
set(gca,'XTick',[])
%axis off
caxis([-2 1])

%axes('position',[.4  .96  .25  .03]);%[.4  .66  .25  .03]);
axes('position',[.25  .94  .5  .025]);%[.4  .66  .25  .03]);
%axes('position',[.592  .66  .1175  .02]);%[.127  .66  .1175  .02]);%.235  .02]);
%imagesc([-eps+similarityRand,-similarityClass2,-similarityClass1,mixtureOdorSimilarity])
imagesc([similarityRand,similarityClass2,similarityClass1]);%,mixtureOdorSimilarity])
axis off
caxis([-2 1]); %caxis([-1 1]); %
text( No/2, 1, 'rand','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
text( No+No/2, 1, '30%','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
text( 2*No+No/2, 1, '70%','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
text( 3*No+No/2, 1, 'mix','fontweight','bold','horizontalalignment','center','color',[1 1 1]);




% piriform -----------
%axes('position',[.7  .7  .25  .25])
axes('position',[.25  .62  .5  .15])
%subplot(3,12,(10:11));
pPlotIdx = 1:100;
piriformRepScaledForFig = piriformRep(pPlotIdx,oPlotIdx)/max(.7*max(piriformRep(pPlotIdx,oPlotIdx)));
imagesc(piriformRepScaledForFig);
caxis([-2 1])


box off;
xlabel('Odor','fontsize',12)
ylabel('Piriform Neuron','fontsize',12); %['Piriform Neuron (',num2str(Ny),' total)']
%title('Piriform','fontsize',12);
set(gca,'YTick',[])
set(gca,'XTick',[])
%axis off
caxis([-2 1])

% axes('position',[.7  .96  .25  .03]); %[.7  .66  .25  .03]);
% %axes('position',[.713  .66  .1175  .02]); %[.398  .66  .1175  .02]);
% %imagesc([-eps+similarityRand,-similarityClass2,-similarityClass1,mixtureOdorSimilarity])
% imagesc([similarityRand,similarityClass2,similarityClass1]);%,mixtureOdorSimilarity])
% axis off
% caxis([-2 1])
% text( No/2, 1, 'rand','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
% text( No+No/2, 1, '30%','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
% text( 2*No+No/2, 1, '70%','fontweight','bold','horizontalalignment','center','color',[1 1 1]);
% %text( 3*No+No/2, 1, 'mix','fontweight','bold','horizontalalignment','center','color',[1 1 1]);








% OVERLAP -----------------------------------------
axes('position',[.55  .33  .35  .2]); hold on
%axes('position',[.55  .08  .37  .37]); hold on
%axes('position',[.55  .1  .4  .4])
%subplot(3,12,[20:23,32:35]);
%[~,randomizer] = sort(rand(size(z.pirRepFreqRandExpected)),'descend');
%for k=1:length(z.pirRepFreqRandExpected)
%    j=randomizer(k);
    plot(z.pirRepFreqClass2Expected,z.pirRepFreqClass2Observed,'.','markersize',7,'color',class2Color);hold on
    plot(z.pirRepFreqClass1Expected,z.pirRepFreqClass1Observed,'.','markersize',7,'color',class1Color);
    plot(z.pirRepFreqRandExpected,z.pirRepFreqRandObserved,'.','markersize',7,'color',[0 0 0]); hold on
%end
%plot(pirRepFreqRandExpected(similarityIsClassIdx),pirRepFreqRandObserved(similarityIsClassIdx),'.','color',[0 .7 0],'markersize',8)
box off;
xlabel('Chance Frequency','fontsize',12)
ylabel('Observed Frequency','fontsize',12)
xlim([0 .03]);
ylim([0 .03]);
set(gca,'XTick',0:.01:.03,'fontsize',12);
set(gca,'YTick',0:.01:.03,'fontsize',12);
plot(xlim,ylim,'color',[.5 .5 .5],'linewidth',1)
alpha(.7); % (transparency)
%saveas(figureUno,'fig1.pdf')




% QUANTIFY SELECTIVITY ---------------------------------------------------
selFrac = 0.6; %how many neurons respond to at least this fraction of class odorants?
class7 = sum(piriformRep(:,2*No+1:end)>0,2) > selFrac*No;
class7sel = sum(class7)/length(class7);
class3 = sum(piriformRep(:,No+1:2*No)>0,2) > selFrac*No;
class3sel = sum(class3)/length(class3);
classRand = sum(piriformRep(:,1:No)>0,2) > selFrac*No;
classRandsel = sum(classRand)/length(classRand);

display([num2str(class7sel),' of pir neurons respond to 60/100 odorants with 70% overlap'])
display([num2str(class3sel),' of pir neurons respond to 60/100 odorants with 30% overlap'])
display([num2str(classRandsel),' of pir neurons respond to 60/100 odorants with random overlap'])


display('figure 1 complete')


