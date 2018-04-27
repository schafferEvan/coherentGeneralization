function x = makeXmagnitudes(xtemp,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig)

x = xtemp.*lognrnd(glomActMu,glomActSig,Nx,No);

while sum(sum(x>glomExpMax))>0
    x(x>glomExpMax)=lognrnd(glomActMu,glomActSig,size(x(x>glomExpMax))); %0;
end
%x = x/glomExpMax;



%% old way (exponential distribution)
% x = xtemp*glomExpMin + xtemp.*random('exp',glomExpMean,Nx,No); %x = x.*random('exp',glomExpMean,Nx,No);
% 
% while sum(sum(x>glomExpMax))>0
%     x(x>glomExpMax)=random('exp',glomExpMean,size(x(x>glomExpMax))); %0;
% end
% %x = x/glomExpMax;