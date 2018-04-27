
figure;
l = {};
d = dir([pwd,'/paramSweep/']);
for i=1:length(d)
    if isempty(strfind(d(i).name,'Sx'))
        continue
    end
    if ~isempty(strfind(d(i).name,'Sx=0.1'))
        continue
    end
    d2 = dir([pwd,'/paramSweep/',d(i).name,'/fly*']);
    l{length(l)+1} = d(i).name;
    
    popTot = load([pwd,'/paramSweep/',d(i).name,'/',d2(1).name]);
    rCCmean = nanmean(popTot(1).rCC,3);
    for j=2:length(d2)
        popTot(j) = load([pwd,'/paramSweep/',d(i).name,'/',d2(j).name]);
        rCCmean = [rCCmean,nanmean(popTot(j).rCC,3)];
    end
    
    dSteps = popTot(1).dSteps; %[popTotS(1).dSteps, popTot(1).dSteps(2:end)]; %[t1.dSteps,t101.dSteps];
    rCCmean = nanmean(rCCmean(:,:),2);
    
    semilogx( dSteps, rCCmean, 'linewidth',1); hold on
    xlabel('Number of Kenyon Cells','fontsize',11)
    ylabel('MBON Correlation','fontsize',11)
    %[hl,legHndl] = legend('SNR','Correlation','Accuracy','Agreement (\theta=0.5)','Agreement (\theta=0.9)','Location','Southeast');
    set(gca,'XTick',[1e1,1e2,1e3,1e4])
    set(gca,'XTickLabel',{'10^1','10^2','10^3','10^4'})
    ylim([0 1])
    xl = xlim;
    xlim([1e2 2*10^4])
    box off
    drawnow
end
legend(l)