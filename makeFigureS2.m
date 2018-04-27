

paperFig2 = figure; 
set(gcf,'color','w')
set(gcf,'Units','Centimeters')
%pos=get(gcf,'Position');
%pos(1)=1000;
pos(3)=11.4; %17.4; %8.5;%pos(3)*1.5;
pos(4)=14.0;%pos(4)*1.5;
set(gcf,'Position',pos)

% annotation('textbox', [.52 .37 .05 .05],...
%     'String', 'F','BackgroundColor','none','Color','k',...
%     'LineStyle','none','fontsize',16,'HorizontalAlignment','Center');
%fiAxes = axes('position',[.61  .12  .13  .19]); 
%fiiAxes = axes('position',[.84  .12  .13  .19]); 
fiAxes = axes('position',[.12  .25  .38  .5]); 
fiiAxes = axes('position',[.55  .25  .38  .5]); 


%% load large values for POP

%l4f = load('/Users/evan/Dropbox/_AxelLab/matlab/_dan/_pop4thOrderResp_0.0/_v2_forPop/pop4thOrderResp_POP_0_iter_73887.mat');

popDirPOP = [pwd,'/_simResults/_pirSims/pop/big/'];
popDLS = dir([popDirPOP,'pop4thOrderResp*']);
popTotPop = load([popDirPOP,popDLS(1).name]);
szp = size(popTotPop(1).rAgr(1,:,:,:,:));
popAgrVec = zeros( szp(1), szp(2), szp(3), szp(4), szp(5)*length(popDLS) );
%popAgrUnVec = zeros( szp(1), szp(2), szp(3), szp(4), szp(5)*length(popTotPop) );
popAgrVec(:,:,:,:,1:szp(5)) = mean(popTotPop(1).rAgr,1); %popTotPop(1).rAgr(1,:,:,:,:);
%popAgrUnVec(:,:,:,:,1:szp(5)) = popTotPop(1).rAgrUn(1,:,:,:,:);
for j=2:length(popDLS)
    popTotPop(j) = load([popDirPOP,popDLS(j).name]);
    popAgrVec(:,:,:,:,(j-1)*szp(5)+(1:szp(5))) = mean(popTotPop(j).rAgr,1); %popTotPop(j).rAgr(1,:,:,:,:);
    %popAgrUnVec(:,:,:,:,(j-1)*szp(5)+(1:szp(5))) = popTotPop(j).rAgrUn(1,:,:,:,:);
end
rCCpMean50 = squeeze(nanmean( popAgrVec, 5));
%if ~isfield(popTotPop(1),'pSteps')
%    popTotPop(1).pSteps = .5:.05:.95;
%end

%% load small values for POP
popDirS = [pwd,'/_simResults/_pirSims/pop/small/'];
popDLS = dir([popDirS,'pop4thOrderResp*']);
popTotPopS = load([popDirS,popDLS(1).name]);
szpS = size(popTotPopS(1).rAgr(1,:,:,:,:));
popAgrVecS = zeros( szpS(1), szpS(2), szpS(3), szpS(4), szpS(5)*length(popDLS) );
%popAgrUnVec = zeros( szp(1), szp(2), szp(3), szp(4), szp(5)*length(popTotPop) );
popAgrVecS(:,:,:,:,1:szpS(5)) = mean(popTotPopS(1).rAgr,1); %popTotPop(1).rAgr(1,:,:,:,:);
%popAgrUnVec(:,:,:,:,1:szp(5)) = popTotPop(1).rAgrUn(1,:,:,:,:);
for j=2:length(popDLS)
    popTotPopS(j) = load([popDirS,popDLS(j).name]);
    popAgrVecS(:,:,:,:,(j-1)*szpS(5)+(1:szpS(5))) = mean(popTotPopS(j).rAgr,1); %popTotPop(j).rAgr(1,:,:,:,:);
    %popAgrUnVec(:,:,:,:,(j-1)*szp(5)+(1:szp(5))) = popTotPop(j).rAgrUn(1,:,:,:,:);
end
rCCpMean50S = squeeze(nanmean( popAgrVecS, 5));
smLgCutPt = 14;

rCCpMean50ag = cat(1, rCCpMean50S, rCCpMean50(smLgCutPt:end,:,:));
dStepsAg = [popTotPopS(1).dSteps,popTotPop(1).dSteps(smLgCutPt:end)];


popAgrNorm = nan(size(rCCpMean50ag));
popAgrNormBase = nan(size(rCCpMean50ag));
for i=1:size(rCCpMean50ag,1)
    for j=1:size(rCCpMean50ag,2)
        for k=1:size(rCCpMean50ag,3)
            % Base is expected value of an untrained population (Binomial):
            Bp=0.5; BN=j; Bx=popTotPop(1).pSteps(k)*BN; 
            popAgrNormBase(i,j,k) = binocdf(BN-Bx,BN,Bp) + binocdf(Bx-1,BN,Bp,'upper'); 
            
            % Normalize agreement by 1-Base, which is the dynamic range:
            popAgrNorm(i,j,k) = (rCCpMean50ag(i,j,k)-popAgrNormBase(i,j,k)) / (1-popAgrNormBase(i,j,k));
        end
    end
end



%annotation('textbox', [.57 .5433 .05 .05],...
axes(fiAxes)
c = .85*parula(size(popAgrNorm(:,end,2:end-1),3));
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');
colormap(c)
semilogx( 10*dStepsAg, squeeze(popAgrNorm(:,10,2:end-1)),'linewidth',1)
hold on; 
semilogx( 2*dStepRef, agrRef,'--','linewidth',1,'color',[.5 .5 .5])
ylim([0 1.2]); xlim([1e4 1e7])
set(gca,'XTick',[1e4,1e5,1e6,1e7])
set(gca,'XTickLabel',{'','10^5','','10^7'})
set(gca,'YTick',0:.25:1)
set(gca,'YTickLabel',{'','','0.5','','1.0'})
box off
xlabel({'Total number of'; 'piriform neurons'},'fontsize',11); 
ylabel('Pop. Agreement','fontsize',11); 
cb = colorbar('Location','north');
cb.Position(4) = 0.02; cb.Position(2) = 0.33; cb.FontSize = 8;
cb.Ticks = [0 1]; cb.TickLabels = {popTotPop(1).pSteps(2),popTotPop(1).pSteps(end-1)};
set(get(cb,'title'),'string','\phi (N_z=10)');

axes(fiiAxes)
c = .85*parula(size(popAgrNorm(:,end,2:end-1),3));
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');
colormap(c)
semilogx( popTotPop(1).N4*dStepsAg, squeeze(popAgrNorm(:,end,2:end-1)),'linewidth',1)
hold on
semilogx( 2*dStepRef, agrRef,'--','linewidth',1,'color',[.5 .5 .5])
ylim([0 1.2]); xlim([1e4 1e7])
set(gca,'XTick',[1e4,1e5,1e6,1e7])
set(gca,'XTickLabel',{'','10^5','','10^7'})
set(gca,'YTick',0:.25:1)
set(gca,'YTickLabel',{'','','0.5','','1.0'})
box off
xlabel({'Total number of'; 'piriform neurons'},'fontsize',11); 
ylabel('Pop. Agreement','fontsize',11); 
cb = colorbar('Location','north');
cb.Position(4) = 0.02; cb.Position(2) = 0.33; cb.FontSize = 8;
cb.Ticks = [0 1]; cb.TickLabels = {popTotPop(1).pSteps(2),popTotPop(1).pSteps(end-1)};
set(get(cb,'title'),'string',['\phi (N_z=',num2str(popTotPop(1).N4),')']);
