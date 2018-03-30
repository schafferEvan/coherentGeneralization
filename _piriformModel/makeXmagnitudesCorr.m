function x = makeXmagnitudesCorr(xtemp,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig,glomCorr)

if length(glomCorr)==1
    S = glomCorr*(ones(No)-eye(No))+eye(No);
else
    S = glomCorr;
end
valTemp = MvLogNRand(glomActMu*ones(1,No),glomActSig*ones(1,No),Nx,S);
%corrcoef(valTemp(:,1),valTemp(:,2))

x = xtemp.*valTemp;

% while sum(sum(x>glomExpMax))>0
%     x(x>glomExpMax)=lognrnd(glomActMu,glomActSig,size(x(x>glomExpMax))); %0;
% end
