
x = xlsread('/Users/evan/Dropbox/_AxelLab/matlab/_dan/_data/_fig1distanceStats/nearest neighbor 180120 compiled.xlsx');
figure; loglog(x(:,1),x(:,2),'k.','markersize',4)
xl = xlim; yl = ylim;
xlim([min(xl(1),yl(1)) max(xl(2),yl(2))])
ylim([min(xl(1),yl(1)) max(xl(2),yl(2))])
xl = xlim; yl = ylim;
hold on; plot([xl(1), xl(2)],[yl(1), yl(2)],'linewidth',1,'color',[.8 .8 .8])
box off
set(gca,'Xtick',[20,40,80]); %0:20:100)
set(gca,'Ytick',[20,40,80]); %0:20:100)
xlabel({'Mean Random Shuffled'; ['N.N. Dist. (\mu','m)']})
ylabel({'Mean Observed'; ['N.N. Dist. (\mu','m)']})
set(gcf,'color','w')
set(gcf,'Units','Centimeters')
pos(4)=4;
pos(3)=4*(xl(2)/yl(2));
set(gcf,'Position',pos)
set(gca,'fontsize',8)


pVal = signrank(x(:,1),x(:,2));
display(['p = ',num2str(pVal)])