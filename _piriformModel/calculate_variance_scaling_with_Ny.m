
Nx = 500;
NyList = 500:500:10000;
c = zeros(size(NyList));
v = zeros(size(NyList));
iters = 100;

parfor j=1:length(NyList)
    display( j/length(NyList) )
    
    tempc = 0;
    tempv = 0;
    for i = 1:10
        xa = randn(Nx,1)/Nx;
        xb = randn(Nx,1)/Nx;
        J1 = randn(NyList(j),Nx)/NyList(j);
        J2 = randn(NyList(j),Nx)/NyList(j);
        tempc = tempc + xa'*(J1'*J1)*xb*xa'*(J2'*J2)*xb;
        tempv = tempv + xa'*(J1'*J1)*xb*xa'*(J1'*J1)*xb;
    end
    
    c(j) = tempc/iters;
    v(j) = tempv/iters;
end

figure; semilogy(NyList,c,'.'); hold all; semilogy(NyList,v,'.')
legend({'correlation across mice','variance within mice'})
    