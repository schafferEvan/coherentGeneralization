
% class and mix replaced by "correlated" on 12/27/15

global figureTwo

fig2sup = figure; 
figureTwo = figure;

set(gcf,'color','w')
pos=get(gcf,'Position');
pos(3)=pos(3);%*.75;
pos(4)=pos(4)*1.75;
set(gcf,'Position',pos)


[~,~,~,~,~,~,plottableCorrSortMono,plottableCorrSortMix,plotComponentIsShared,ex1,ex2,odorLabels] = plotDataCorrelation;

% ------------------------------------------------
figure(figureTwo);

%axes('position',[.2  .8  .2  .2*pos(3)/pos(4)]); hold all
axes('position',[.25  .78  .2  .2*pos(3)/pos(4)]); hold all
valid1 = (ex1.respA>0) & (ex1.respB>0);
plot(ex1.respA(valid1),ex1.respB(valid1),'ko','markersize',4)
title('Hex vs Oct','fontsize',12);
xlabel('Hexanal Resp. (\Delta F/F)','fontsize',12);
ylabel('Octanal Resp. (\Delta F/F)','fontsize',12);
%annotation('textbox', [.32 .9 .1 .02],...
%    'String', ['CC = ',num2str(1/100*round(100*ex1.cc))],'BackgroundColor','none','Color','k',...
%    'LineStyle','none','fontsize',10,'HorizontalAlignment','Center');
text(.25, .05, ['CC = ',num2str(1/100*round(100*ex1.cc))],'BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',10,'HorizontalAlignment','Center');
annotation('textbox', [.15 .94 .05 .05],...
    'String', 'A','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',32,'HorizontalAlignment','Center');


%axes('position',[.6  .8  .2  .2*pos(3)/pos(4)]); hold all
axes('position',[.65  .78  .2  .2*pos(3)/pos(4)]); hold all
valid2 = (ex2.respA>0) & (ex2.respB>0);
plot(ex2.respA(valid2),ex2.respB(valid2),'ko','markersize',4)
title('EMS+Oct vs Oct','fontsize',12);
xlabel('EMS+Oct Resp. (\Delta F/F)','fontsize',12);
ylabel('Octanal Resp. (\Delta F/F)','fontsize',12);
%annotation('textbox', [.72 .9 .1 .02],...
%    'String', ['CC = ',num2str(1/100*round(100*ex2.cc))],'BackgroundColor','none','Color','k',...
%    'LineStyle','none','fontsize',10,'HorizontalAlignment','Center');
text(.45, .1, ['CC = ',num2str(1/100*round(100*ex2.cc))],'BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',10,'HorizontalAlignment','Center');
annotation('textbox', [.55 .94 .05 .05],...
    'String', 'B','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',32,'HorizontalAlignment','Center');


% 
% M = max(dataF0(bothActiveF0,1));
% m = min(dataF0(bothActiveF0,1));
% dataF0(:,1) = (dataF0(:,1)-m) * (numBins-1) / (M-m) + 1;
% M = max(dataF0(bothActiveF0,2));
% m = min(dataF0(bothActiveF0,2));
% dataF0(:,2) = (dataF0(:,2)-m) * (numBins-1) / (M-m) + 1;
% plot(dataF0(bothActiveF0,1),dataF0(bothActiveF0,2),'w.','markersize',4)
% set(gca,'XTick',[]); set(gca,'YTick',[])
% ylabel('Odor B','fontsize',12);%,'FontWeight','bold'); 
% title('Example 1','fontsize',12);
% xlim([0 numBins])
% ylim([0 numBins])


% ------------------------------------------------

%plotModelCorrVsInput;
figure(figureTwo);

%classStruct = load('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatClass.mat');%,'corrMat','f');
%mixStruct = load('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatMix.mat');%,'corrMat','f');
corrStruct = load('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatCorr.mat');%,'corrMat','f');

% % rejection plot v1 (background color)
% axes('position',[.25  .1  .65  .2]); hold on
% %imagesc(corrStruct.fracSignificant')
% imagesc(corrStruct.fracSignificant')
% colormap(bone); brighten(.2)
% %colorbar;
% axis off

% rejection plot v2 (area plot)
axes('position',[.1  .05  .15  .2]); hold on
fill([0,-1+corrStruct.fracSignificant,0],[corrStruct.f(1),corrStruct.f,corrStruct.f(end)],[.92 .92 .92],'linewidth',1); hold on
xlim([-1 0])
ylim([0 1])
%axis off
set(gca,'Xtick',{})
set(gca,'Ytick',{})
set(gca,'YAxisLocation','right')
box off;
xlabel('P(reject)','fontsize',12)
%text(-.8,-.15,'P(reject)','fontsize',12)




%axes('position',[.15  .1  .65  .2]); hold on
axes('position',[.25  .05  .65  .2]); hold on
set(gca,'XDir','reverse')

fill( [corrStruct.f,fliplr(corrStruct.f)], [mean(corrStruct.corrMat)-std(corrStruct.corrMat), fliplr(mean(corrStruct.corrMat)+std(corrStruct.corrMat))] , 'g','edgecolor','none'); alpha(.15);
plot( corrStruct.f, mean(corrStruct.corrMat) , 'g','linewidth',1 )
%errorbar( corrStruct.f, mean(corrStruct.corrMat) ,std(corrStruct.corrMat) , 'g.' )
% %errorbar(classStruct.f, mean(classStruct.corrMat) , std(classStruct.corrMat) , 'g.' )
errorbar(corrStruct.f(3), mean(corrStruct.corrMat(:,3)) , std(corrStruct.corrMat(:,3)) , 'k.' );%, 'linewidth',1.5)

%legend('Mix','Class','Random', 'Location','Southeast')
xlabel('fraction input overlap','fontsize',12);
ylabel('Model Correlation Coefficient (mean \pm \sigma)','fontsize',12);
ylim([-.2 1])
xlim([0 1])
box off
view(90,90)
set(gca,'color','none')
annotation('textbox', [.15 .25 .05 .05],...
    'String', 'D','BackgroundColor','none','Color','k',...
    'LineStyle','none','fontsize',32,'HorizontalAlignment','Center');


% figure; hist(corrStruct.corrMat(:,3),20)
% xlim([-.2 1])
% xlabel('Correlation Coefficient','fontsize',12);
% ylabel('Probability Density','fontsize',12);
% 
% figure(fig2sup)
% isGreaterThanMean = plottableCorrSortMono>mean(corrStruct.corrMat(:,3));
% monoRandpVal = zeros(size(plottableCorrSortMono));
% for j=1:length(monoRandpVal)
%     if isGreaterThanMean(j)
%         monoRandpVal(j) = sum( corrStruct.corrMat(:,3)>=plottableCorrSortMono(j) ) / length( corrStruct.corrMat(:,3) );
%     else
%         monoRandpVal(j) = sum( corrStruct.corrMat(:,3)<=plottableCorrSortMono(j) ) / length( corrStruct.corrMat(:,3) );
%     end
% end
% markSize = 5;
% stem(1:length(plottableCorrSortMono),monoRandpVal,'k','markersize',markSize); 
% set(gca,'XTick',1:length(odorLabels),'XTickLabel',odorLabels);%,'XTickLabelRotation',-90);
% ylabel('One-sided p-value','fontsize',12)
% view(90,90)
% box off
% 
% szCmat = size(corrStruct.corrMat);
% modelRandMean = mean(corrStruct.corrMat(:,3));
% modelRandStd = std(corrStruct.corrMat(:,3));
% modelClassMean = mean(corrStruct.corrMat);
% modelClassStd = std(corrStruct.corrMat);
% modelMixMean = mean(corrStruct.corrMat);
% modelMixStd = std(corrStruct.corrMat);
% odorVals = [plottableCorrSortMono;plottableCorrSortMix];
% 
% modelRandStds = abs(odorVals-modelRandMean)/modelRandStd;
% modelClassStds = NaN*ones(length(odorVals),szCmat(2));
% modelMixStds = NaN*ones(length(odorVals),szCmat(2));
% for j=1:szCmat(2)
%     modelClassStds(:,j) = abs(odorVals-modelClassMean(j))/modelClassStd(j);
%     modelMixStds(:,j) = abs(odorVals-modelMixMean(j))/modelMixStd(j);
% end


% % -- OLD PANEL "E" -------------
% figure(figureTwo);
% axes('position',[.83  .35  .12  .4]); hold all
% stem(1:length(odorVals),modelRandStds,'k','markersize',markSize); 
% stem(1:length(odorVals),min(modelClassStds,[],2),'g','markersize',markSize); 
% stem(1:length(odorVals),min(modelMixStds,[],2),'b','markersize',markSize); 
% set(gca,'XTick',1:length(odorLabels),'XTickLabel',[]);%,'XTickLabelRotation',-90);
% %set(gca,'XTick',1:length(odorLabels),'XTickLabel',odorLabels);%,'XTickLabelRotation',-90);
% ylabel('\sigma','fontsize',12)
% %ylabel('Standard deviations from the mean','fontsize',12)
% view(90,90)
% box off

