
function calculateImportanceFig_pInv
global Nx Ny Sx Sc glomExpMean glomExpMax Sp NoR NoC NoM classOdorOverlap N4 S4 S4c templateOdor paperFig




% *************************************************************************
% *************************************************************************
% *************************************************************************

%              DEPRECATED.  REPLACED BY
%              "calculateImportanceFigParts_pInv.m"

% *************************************************************************
% *************************************************************************
% *************************************************************************









NoR6 = 5*NoR;%*50; % NOTE: NoR6 and NoC6 must be the same
NoC6 = 5*NoC;%*50;
NoM6 = 100;
numOdorsToPlot = 200;
normalizeResp = false;

% ODORS ---------------------------------------------------------

[glomRepTempRand,odorSimilarity] = makeOdors('classWithoutTemplate',0,NoR6);
similarityMeasureFull = -eps*(odorSimilarity<0);
[glomRepTempClass,odorSimilarity] = makeOdors('classWithTemplate',1,NoC6); 
similarityMeasureFull = [-2+eps*(odorSimilarity<0),similarityMeasureFull];
[glomRepTempMix,odorSimilarity] = makeOdors('mix',0,NoM6);
mixtureOdorSimilarity = odorSimilarity;
similarityMeasureFull = [similarityMeasureFull, -1+odorSimilarity/NoM6];
glomerulusRep = [glomRepTempClass, glomRepTempRand, glomRepTempMix];

glomRepMouse2 = makeOdors('classWithoutTemplate',0,NoR6);

% PIRIFORM -------------------------------------------------------

[piriformRep,~,~,J,maxy] = makePiriform(glomerulusRep);  % mouse 1
pirClassRand = makePiriform(glomRepMouse2,Sp,J,maxy); % mouse 2
%pirClassRand = pirClassRand(:,NoC6+1:NoC6+NoR6);


%%  CLASS DISCRIMINANT & IMPORTANCE (class example) -------------------------

pirClass = piriformRep(:,1:NoC6);
pirRand = piriformRep(:,NoC6+1:NoC6+NoR6);
mc = mean(pirClass,2);
mr = mean(pirRand,2);
%mcMat = repmat( mc , 1,200);
%mrMat = repmat( mr , 1,200);
%withinCov = (pirClass-mcMat)*(pirClass'-mcMat') + (pirRand-mrMat)*(pirRand'-mrMat');
%FisherLDweights = withinCov \ (mc-mr);
%FisherLDweights = FisherLDweights/abs(sum(FisherLDweights));

importance = (mc-mr);% ./ mr;% ./ (mc+mr);

[~,importanceLoc] = sort(importance,'descend');
importanceBound1 = 40000;%500;
importanceBound2 = find(importance(importanceLoc)>0, 1, 'last'); %Ny;
importanceSparseBound = 0.1; %.99;
temp = rand(size(mc))<importanceSparseBound;

A = full([pirRand,pirClass]);
mtemp = mean(mean(A));
classWeightsFull = [-ones(1,NoR6),ones(1,NoC6)]*pinv(A-mtemp);
%classWeightsFull = (mc-mr).*( (mc-mr)>0 );        
classWeightsNoTop = classWeightsFull;
classWeightsNoTop(importanceLoc(1:importanceBound1)) = 0;
classWeightsNoBottom = classWeightsFull;
classWeightsNoBottom(importanceLoc(importanceBound1+1:importanceBound2)) = 0;
classWeightsRand = classWeightsFull;
classWeightsRand(temp) = 0;


classFourthOrderFull = classWeightsFull*piriformRep;
classFourthOrderNoTop = classWeightsNoTop*piriformRep;
classFourthOrderNoBottom = classWeightsNoBottom*piriformRep;
classFourthOrderRand = classWeightsRand*piriformRep;

if normalizeResp
    %normalization & shifting
    classFourthOrderFull = (classFourthOrderFull-min(classFourthOrderFull))/(max(classFourthOrderFull)-min(classFourthOrderFull));
    classFourthOrderNoTop = (classFourthOrderNoTop-min(classFourthOrderNoTop))/(max(classFourthOrderNoTop)-min(classFourthOrderNoTop));
    classFourthOrderNoBottom = (classFourthOrderNoBottom-min(classFourthOrderNoBottom))/(max(classFourthOrderNoBottom)-min(classFourthOrderNoBottom));
    classFourthOrderRand = (classFourthOrderRand-min(classFourthOrderRand))/(max(classFourthOrderRand)-min(classFourthOrderRand));
else
    %just shifting
    classFourthOrderFull = classFourthOrderFull-min(classFourthOrderFull);
    classFourthOrderNoTop = classFourthOrderNoTop-min(classFourthOrderNoTop);
    classFourthOrderNoBottom = classFourthOrderNoBottom-min(classFourthOrderNoBottom);
    classFourthOrderRand = classFourthOrderRand-min(classFourthOrderRand);
end

%% CLASS DISCRIMINANT & IMPORTANCE (RANDOM example) -------------------------

%pirClassRand = pirClassRand;%(:,NoC6+1:NoC6+NoR6);
%pirRand = piriformRep(:,NoC6+1:NoC6+NoR6);
mcr = mean(pirClassRand,2);
%mr = mean(pirRand,2);

importanceRand = (mcr-mr);% ./ mr;% ./ (mcr+mr);

[~,importanceLocRand] = sort(importanceRand,'descend');
importanceBound1Rand = 40000;%500;
importanceBound2Rand = find(importanceRand(importanceLocRand)>0, 1, 'last'); %Ny;
%importanceSparseBound = Sc; %.99;
temp = rand(size(mcr))<importanceSparseBound;

A = full([pirRand,pirClassRand]);
mtemp = mean(mean(A));
classWeightsFullRand = [-ones(1,NoR6),ones(1,NoR6)]*pinv(A-mtemp);
%classWeightsFullRand = (mcr-mr).*( (mcr-mr)>0 );           
classWeightsNoTopRand = classWeightsFullRand;
classWeightsNoTopRand(importanceLocRand(1:importanceBound1Rand)) = 0;
classWeightsNoBottomRand = classWeightsFullRand;
classWeightsNoBottomRand(importanceLocRand(importanceBound1Rand+1:importanceBound2Rand)) = 0;
classWeightsRandRand = classWeightsFullRand;
classWeightsRandRand(temp) = 0;

pirRepTemp = piriformRep;
pirRepTemp(:,1:NoC6) = pirClassRand;%(:,NoC6+1:NoC6+NoR6);
classFourthOrderFullRand = classWeightsFullRand*pirRepTemp;
classFourthOrderNoTopRand = classWeightsNoTopRand*pirRepTemp;
classFourthOrderNoBottomRand = classWeightsNoBottomRand*pirRepTemp;
classFourthOrderRandRand = classWeightsRandRand*pirRepTemp;

if normalizeResp
    %normalization & shifting
    classFourthOrderFullRand = (classFourthOrderFullRand-min(classFourthOrderFullRand))/(max(classFourthOrderFullRand)-min(classFourthOrderFullRand));
    classFourthOrderNoTopRand = (classFourthOrderNoTopRand-min(classFourthOrderNoTopRand))/(max(classFourthOrderNoTopRand)-min(classFourthOrderNoTopRand));
    classFourthOrderNoBottomRand = (classFourthOrderNoBottomRand-min(classFourthOrderNoBottomRand))/(max(classFourthOrderNoBottomRand)-min(classFourthOrderNoBottomRand));
    classFourthOrderRandRand = (classFourthOrderRandRand-min(classFourthOrderRandRand))/(max(classFourthOrderRandRand)-min(classFourthOrderRandRand));
else
    %just shifting
    classFourthOrderFullRand = classFourthOrderFullRand-min(classFourthOrderFullRand);
    classFourthOrderNoTopRand = classFourthOrderNoTopRand-min(classFourthOrderNoTopRand);
    classFourthOrderNoBottomRand = classFourthOrderNoBottomRand-min(classFourthOrderNoBottomRand);
    classFourthOrderRandRand = classFourthOrderRandRand-min(classFourthOrderRandRand);
end



%% IMPORTANCE FIGURE
figure; 
set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=pos(3); %*3/4;
%pos(4)=pos(4)*2;
pos(4)=pos(4)*1;
set(gcf,'Position',pos)
%respBins = 0:.05:1;
color1 = [1 .6 .3];
color2 = [1 1 .4];
classColor = [0 .8 0];
colormap(bone)
maxBins = 1;
%ha = tight_subplot(5, 4, .05, .05, .05); 




%%  PART 1: SENSIBLE CLASS DEFINITIONS ------------

% annotation('textbox', [0 .5 .1 .1],...
%            'String', 'Piriform','BackgroundColor',[0 0 0],'Color',[1 1 1],...
%            'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');
% annotation('textbox', [0 .15 .1 .1],...
%            'String', '4th Order','BackgroundColor',[0 0 0],'Color',[1 1 1],...
%            'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');


%axes(ha(2)); hold on
%axes('position',[.2375  .75  .15  .15*pos(3)/pos(4)]); hold on
axes('position',[.2625  .75  .2  .15*pos(3)/pos(4)]); hold on
%subplot(5,4,2); hold on
odorspace = zeros(20);
odorspace(3:8,12:17) = 1;
imagesc(odorspace)
axis off

annotation('textbox', [.1875 .87 .3 .1],...
    'String', 'Class = Similar Odors','BackgroundColor','none','Color','k',...
    'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');


% axes('position',[.1  .61  .3750  .12]); hold on
% %subplot(5,4,4+(1:2)); hold on
% %plot(importance(importanceLoc(1:1.5*Ny*Sp)),'k','linewidth',1.5)
% x1=[1:importanceBound1,importanceBound1:-1:1];
% x2=[importanceBound1:importanceBound2,importanceBound2:-1:importanceBound1];
% x3=[importanceBound2:Ny,Ny:-1:importanceBound2];
% area(x1,[importance(importanceLoc(1:importanceBound1));0*importance(importanceLoc(1:importanceBound1))],'Facecolor',color1)
% area(x2,[importance(importanceLoc(importanceBound1:importanceBound2));0*importance(importanceLoc(importanceBound1:importanceBound2))],'Facecolor',color2)
% area(x3,[importance(importanceLoc(importanceBound2:Ny));0*importance(importanceLoc(importanceBound2:Ny))],'Facecolor',[.95,.95,.95])
% box off
% %xlim([0 1.5*Ny*Sp]);
% xlim([1 Ny]);
% %ylim([0 1]);
% xlabel('Piriform Neuron Rank','fontsize',10)
% ylabel('Class Selectivity','fontsize',10)

%axes('position',[.1  .44  .3750  .12]); hold on
axes('position',[.15  .55  .3750  .12]); hold on
%subplot(5,4,2*4+[1 2]); hold on
selectedNeuron = importanceLoc(importanceBound1);
bar(numOdorsToPlot+1:2*numOdorsToPlot,pirRand(selectedNeuron,1:numOdorsToPlot),'k','EdgeColor','none');
bar(1:numOdorsToPlot,pirClass(selectedNeuron,1:numOdorsToPlot),'FaceColor', classColor,'EdgeColor','none');
box off
ylim([0 .25]);
xlim([1, 2*numOdorsToPlot+1]);
xlabel('Odor Identity','fontsize',10)
ylabel('Response Magnitude','fontsize',10)
legend('Random','Class','Location',[1.5 .5 .3 .3]);
text(100,.19,num2str( sum( pirClass(selectedNeuron,1:numOdorsToPlot)>0)/numOdorsToPlot ), 'fontweight','bold','color',classColor);
text(300,.19,num2str( sum( pirRand(selectedNeuron,1:numOdorsToPlot)>0)/numOdorsToPlot ), 'fontweight','bold','color','k');
text(-100,.125,'Piriform','Rotation',90,'Backgroundcolor',[.4 .4 .4],'Color',[1 1 1],...
    'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');
text(-100,-0.6,'4th Order','Rotation',90,'Backgroundcolor',[.4 .4 .4],'Color',[1 1 1],...
    'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');

       

axes('position',[.15  .295  .1625  .12]); hold on
%subplot(5,4,4+9); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderFull(NoC6+1:end-NoM6)/classFourthOrderFull(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderFull(1:NoC6)/classFourthOrderFull(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderFull(NoC6+1:end-NoM6)/max(classFourthOrderFull),respBins)/NoR6;
    hc = hist(classFourthOrderFull(1:NoC6)/max(classFourthOrderFull),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderFull(NoC6+1:end-NoM6)), max(classFourthOrderFull(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderFull(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderFull(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
ylabel('Probability','fontsize',10)
%xlabel('Response','fontsize',10)
title('Full','fontsize',12)
%legend('Random','Class');

axes('position',[.3625  .295  .1625  .12]); hold on
%subplot(5,4,4+10); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderRand(NoC6+1:end-NoM6)/classFourthOrderRand(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderRand(1:NoC6)/classFourthOrderRand(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderRand(NoC6+1:end-NoM6)/max(classFourthOrderRand),respBins)/NoR6;
    hc = hist(classFourthOrderRand(1:NoC6)/max(classFourthOrderRand),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderRand(NoC6+1:end-NoM6)), max(classFourthOrderRand(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderRand(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderRand(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
%xlabel('Response','fontsize',10)
title(['Random ',num2str(importanceSparseBound*100),'%'],'fontsize',12)

axes('position',[.15  .1  .1625  .12]); hold on
%subplot(5,4,4+13); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderNoBottom(NoC6+1:end-NoM6)/classFourthOrderNoBottom(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderNoBottom(1:NoC6)/classFourthOrderNoBottom(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderNoBottom(NoC6+1:end-NoM6)/max(classFourthOrderNoBottom),respBins)/NoR6;
    hc = hist(classFourthOrderNoBottom(1:NoC6)/max(classFourthOrderNoBottom),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderNoBottom(NoC6+1:end-NoM6)), max(classFourthOrderNoBottom(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderNoBottom(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderNoBottom(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
ylabel('Probability','fontsize',10)
xlabel('Response','fontsize',10)
title(['Top ',num2str(round(importanceBound1/Ny*100)),'%'],'fontsize',12)
%title(['All But Top ',num2str(importanceBound1),' Ablated'],'fontsize',12)
%set(gca,'color',color1)

axes('position',[.36  .1  .1625  .12]); hold on
%subplot(5,4,4+14); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderNoTop(NoC6+1:end-NoM6)/classFourthOrderNoTop(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderNoTop(1:NoC6)/classFourthOrderNoTop(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderNoTop(NoC6+1:end-NoM6)/max(classFourthOrderNoTop),respBins)/NoR6;
    hc = hist(classFourthOrderNoTop(1:NoC6)/max(classFourthOrderNoTop),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderNoTop(NoC6+1:end-NoM6)), max(classFourthOrderNoTop(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderNoTop(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderNoTop(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
xlabel('Response','fontsize',10)
title(['Bottom ',num2str(round((1-importanceBound1/Ny)*100)),'%'],'fontsize',12)
%title(['Top ',num2str(importanceBound1),' Ablated'],'fontsize',12)
%set(gca,'color',color2)








%%  PART 2: RANDOM CLASS DEFINITIONS ------------


axes('position',[.6875  .75  .2  .15*pos(3)/pos(4)]); hold on
%axes('position',[.7  .5833  1/6  1/6]); hold on
%subplot(5,4,3); hold on
odorspace = zeros(20);
odorlocIdx = 1:20^2;
[~,randOdorLoc] = sort(randn(size(odorlocIdx)),'descend');
odorspace(randOdorLoc(1:20)) = 1;
imagesc(odorspace)
axis off
%title({'Class ='; 'Random Odors'},'fontsize',10)
annotation('textbox', [.61 .87 .35 .1],...
           'String', 'Class = Random Odors','BackgroundColor','none','Color','k',...
           'FontWeight','bold','LineStyle','none','fontsize',14,'HorizontalAlignment','Center');


% axes('position',[.525  .61  .375  .12]); hold on
% %subplot(5,4,6+(1:2)); hold on
% %plot(importanceRand(importanceLocRand(1:1.5*Ny*Sp)),'k','linewidth',1.5)
% x1=[1:importanceBound1Rand,importanceBound1Rand:-1:1];
% x2=[importanceBound1Rand:importanceBound2Rand,importanceBound2Rand:-1:importanceBound1Rand];
% x3=[importanceBound2Rand:Ny,Ny:-1:importanceBound2Rand];
% area(x1,[importanceRand(importanceLocRand(1:importanceBound1Rand));0*importanceRand(importanceLocRand(1:importanceBound1Rand))],'Facecolor',color1)
% area(x2,[importanceRand(importanceLocRand(importanceBound1Rand:importanceBound2Rand));0*importanceRand(importanceLocRand(importanceBound1Rand:importanceBound2Rand))],'Facecolor',color2)
% area(x3,[importanceRand(importanceLocRand(importanceBound2Rand:Ny));0*importanceRand(importanceLocRand(importanceBound2Rand:Ny))],'Facecolor',[.95,.95,.95])
% box off
% %xlim([0 1.5*Ny*Sp]);
% xlim([1 Ny]);
% %ylim([0 1]);
% xlabel('Piriform Neuron Rank','fontsize',10)
% %ylabel('Class Selectivity','fontsize',10)

%axes('position',[.525  .44  .375  .12]); hold on
axes('position',[.575  .55  .375  .12]); hold on
%subplot(5,4,6+(5:6)); hold on
selectedNeuronRand = importanceLocRand(importanceBound1Rand);
bar(numOdorsToPlot+1:2*numOdorsToPlot,pirRand(selectedNeuronRand,1:numOdorsToPlot),'k','EdgeColor','none');
bar(1:numOdorsToPlot,pirClassRand(selectedNeuronRand,1:numOdorsToPlot),'FaceColor', classColor,'EdgeColor','none');
box off
ylim([0 .25]);
xlim([1, 2*numOdorsToPlot+1]);
xlabel('Odor Identity','fontsize',10)
text(100,.19,num2str( sum( pirClassRand(selectedNeuronRand,1:numOdorsToPlot)>0)/numOdorsToPlot ), 'fontweight','bold','color',classColor);
text(300,.19,num2str( sum( pirRand(selectedNeuronRand,1:numOdorsToPlot)>0)/numOdorsToPlot ), 'fontweight','bold','color','k');
%ylabel('Response Magnitude','fontsize',10)
%legend('Random','Class');

axes('position',[.575  .295  .1625  .12]); hold on
%subplot(5,4,6+9); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderFull(NoC6+1:end-NoM6)/classFourthOrderFull(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderFull(1:NoC6)/classFourthOrderFull(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderFullRand(NoC6+1:end-NoM6)/max(classFourthOrderFullRand),respBins)/NoR6;
    hc = hist(classFourthOrderFullRand(1:NoC6)/max(classFourthOrderFullRand),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderFullRand(NoC6+1:end-NoM6)), max(classFourthOrderFullRand(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderFullRand(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderFullRand(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
%xlabel('Response','fontsize',10)
title('Full','fontsize',12)
%legend('Random','Class');

axes('position',[.7875  .295  .1625  .12]); hold on
%subplot(5,4,6+10); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderRand(NoC6+1:end-NoM6)/classFourthOrderRand(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderRand(1:NoC6)/classFourthOrderRand(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderRandRand(NoC6+1:end-NoM6)/max(classFourthOrderRandRand),respBins)/NoR6;
    hc = hist(classFourthOrderRandRand(1:NoC6)/max(classFourthOrderRandRand),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderRandRand(NoC6+1:end-NoM6)), max(classFourthOrderRandRand(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderRandRand(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderRandRand(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
%xlabel('Response','fontsize',10)
title(['Random ',num2str(importanceSparseBound*100),'%'],'fontsize',12)


axes('position',[.575  .1  .1625  .12]); hold on
%subplot(5,4,6+13); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderNoBottom(NoC6+1:end-NoM6)/classFourthOrderNoBottom(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderNoBottom(1:NoC6)/classFourthOrderNoBottom(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderNoBottomRand(NoC6+1:end-NoM6)/max(classFourthOrderNoBottomRand),respBins)/NoR6;
    hc = hist(classFourthOrderNoBottomRand(1:NoC6)/max(classFourthOrderNoBottomRand),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderNoBottomRand(NoC6+1:end-NoM6)), max(classFourthOrderNoBottomRand(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderNoBottomRand(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderNoBottomRand(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
xlabel('Response','fontsize',10)
title(['Top ',num2str(round(importanceBound1Rand/Ny*100)),'%'],'fontsize',12)
%title(['All But Top ',num2str(importanceBound1),' Ablated'],'fontsize',12)
%set(gca,'color',color1)


axes('position',[.7875  .1  .1625  .12]); hold on
%subplot(5,4,6+14); hold on
%plot(glomCorrOneSorted(NoC6+1:end-NoM6),classFourthOrderNoTop(NoC6+1:end-NoM6)/classFourthOrderNoTop(1),'k.','markersize',8)
%plot(glomCorrOneSorted(1:NoC6),classFourthOrderNoTop(1:NoC6)/classFourthOrderNoTop(1),'.', 'color', [0 .8 0],'markersize',8)
if normalizeResp
    respBins = 0:.05:1;
    hr = hist(classFourthOrderNoTopRand(NoC6+1:end-NoM6)/max(classFourthOrderNoTopRand),respBins)/NoR6;
    hc = hist(classFourthOrderNoTopRand(1:NoC6)/max(classFourthOrderNoTopRand),respBins)/NoC6;
else
    maxBins = max( max(classFourthOrderNoTopRand(NoC6+1:end-NoM6)), max(classFourthOrderNoTopRand(1:NoC6)) );
    respBins = 0:maxBins/20:maxBins;
    hr = hist(classFourthOrderNoTopRand(NoC6+1:end-NoM6),respBins)/NoR6;
    hc = hist(classFourthOrderNoTopRand(1:NoC6),respBins)/NoC6;
end
bar(respBins,hr,'k'); hold on
bar(respBins,hc,'FaceColor', classColor)
box off
xlim([0 maxBins]);
ylim([0 0.6]);
%ylabel('Probability','fontsize',10)
xlabel('Response','fontsize',10)
title(['Bottom ',num2str(round((1-importanceBound1Rand/Ny)*100)),'%'],'fontsize',12)
%title(['Top ',num2str(importanceBound1),' Ablated'],'fontsize',12)
%set(gca,'color',color2)



