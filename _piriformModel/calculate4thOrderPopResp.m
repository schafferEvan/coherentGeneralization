
global Ny

% for iterIter = 1:4; batchJob = batch(@calculate4thOrderPopResp,0);
% wait(batchJob); end % doesn't work.  Needs input
%popRespIter = 0;
popRespIter = 1; %popRespIter + 1;

initializeGlobals_Server;
No = 1000; %200;
dSteps = round(logspace(3,6,28)); %[1e4,1e5,1e6]; round(logspace(3,6,10)); %
kSteps = 1;%.1:.2:.9;  %no longer in use
N4 = 200;

tmpTime = clock;
rng(tmpTime(end)); %seed random number generator

glomRepLarge = makeOdors('correlated',0.7,No);
gCC = corrcoef(glomRepLarge);
rCC = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rCCpSum = zeros(length(dSteps), length(kSteps), No, N4); %zeros(Ny/dStep,dimIters);
rCCSingle = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
r = zeros(length(dSteps), length(kSteps), No-1, N4); %zeros(Ny/dStep,dimIters);

rAgr = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rAgrSum = zeros(length(dSteps), length(kSteps), No, N4); %zeros(Ny/dStep,dimIters);
rAgrSingle = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);

REx = zeros(length(dSteps), No-1);
REx2 = zeros(length(dSteps), No-1);
RExSingle = zeros(length(dSteps), No-1);
RExSingle2 = zeros(length(dSteps), No-1);

gCCn = gCC;
for o=1:No
    gCCn((1:No)~=o, o) = gCCn((1:No)~=o, o) - mean(gCCn((1:No)~=o, o));
    gCCn((1:No)~=o, o) = gCCn((1:No)~=o, o) / norm(gCCn((1:No)~=o, o));
end

for j=1:length(dSteps)
    display(j/length(dSteps))
    
    Ny = dSteps(j);
    piriformRepLarge = makePiriform(glomRepLarge);
    piriformRepLarge2 = makePiriform(glomRepLarge);
    
    for k=1:length(kSteps)
        
        % Randomly select sparse connectivity for each 4th order neuron
        WIdx = zeros(N4,Ny);
        for l=1:N4
            yOrder = randperm(Ny);
            %WIdx(l, yOrder( 1:round(Ny*kSteps(k)) ) ) = 1;
            WIdx(l, yOrder( 1:1000 ) ) = 1;
        end
        WIdx2 = zeros(N4,Ny);
        for l=1:N4
            yOrder = randperm(Ny);
            %WIdx(l, yOrder( 1:round(Ny*kSteps(k)) ) ) = 1;
            WIdx2(l, yOrder( 1:1000 ) ) = 1;
        end
        
        for o=1:No
            
            % "train" all 4th order neurons to same odor
            
            W = WIdx .* repmat(piriformRepLarge(:,o)', N4, 1);
            W2 = WIdx2 .* repmat(piriformRepLarge2(:,o)', N4, 1);
            %W = W-mean(mean(W));%(W~=0)));
            W(W~=0) = W(W~=0) - mean(mean(W(W~=0)));
            W2(W2~=0) = W2(W2~=0) - mean(mean(W2(W2~=0)));
            rv = full( W*piriformRepLarge(:,(1:No)~=o) );
            rv2 = full( W2*piriformRepLarge2(:,(1:No)~=o) );
            %rv = rv./repmat(sum(rv,2),1,99);
            R = sum(rv);
            R2 = sum(rv2);
            tmp = R-mean(R);
            tmp = tmp/norm(tmp);
            tmp2 = R2-mean(R2);
            tmp2 = tmp2/norm(tmp2);
            %rCC(j,k,o) = tmp*gCCn((1:No)~=o, o); %tmp(2);
            
            REx(j,:) = tmp;
            REx2(j,:) = tmp2;
            mTh1 = full(median(tmp));
            mTh2 = full(median(tmp2));
            %mTh2 = full(median(gCCn((1:No)~=o, o)));
            %rAgr(j,k,o) = ( sum( (gCCn((1:No)~=o, o)>mTh2) & (tmp>mTh1)'  ) ...
            %    +  sum( (gCCn((1:No)~=o, o)<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
            rAgr(j,k,o) = ( sum( (tmp2>mTh2) & (tmp>mTh1)'  ) ...
                +  sum( (tmp2<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
            
            Nsingle = 10;
            tmpV = zeros(1,floor(N4/Nsingle));
            tmpV2 = zeros(1,floor(N4/Nsingle));
            for l=1:length(tmpV)
                tmp = sum(rv( 1+(l-1)*Nsingle:l*Nsingle, :));
                tmp2 = sum(rv2( 1+(l-1)*Nsingle:l*Nsingle, :));
                tmp = tmp-mean(tmp);
                tmp = tmp/norm(tmp);
                tmp2 = tmp2-mean(tmp2);
                tmp2 = tmp2/norm(tmp2);
                %tmpV(l) = tmp*gCCn((1:No)~=o, o); %tmp(2);
                
                mTh1 = full(median(tmp));
                mTh2 = full(median(tmp2));
                %mTh2 = full(median(gCCn((1:No)~=o, o)));
                %tmpV2 = ( sum( (gCCn((1:No)~=o, o)>mTh2) & (tmp>mTh1)'  ) ...
                %    +  sum( (gCCn((1:No)~=o, o)<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
                tmpV2 = ( sum( (tmp2>mTh2) & (tmp>mTh1)'  ) ...
                    +  sum( (tmp2<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
            end
            %rCCSingle(j,k,o) = mean(tmpV);
            rAgrSingle(j,k,o) = mean(tmpV2);
            RExSingle(j,:) = tmp;
            
            %r(j,k,:,:) = rv'; %this is wrong
            
            Rc = cumsum(rv);
            Rc2 = cumsum(rv2);
            for l=1:N4
                tmp = Rc(l,:)-mean(Rc(l,:));
                tmp = tmp/norm(tmp);
                tmp2 = Rc2(l,:)-mean(Rc2(l,:));
                tmp2 = tmp2/norm(tmp2);
                %tmp = corrcoef( Rc(l,:), gCC((1:No)~=o, o) );
                %rCCpSum(j,k,o,l) = tmp*gCCn((1:No)~=o, o); %tmp(2);
                
                mTh1 = full(median(tmp));
                mTh2 = full(median(tmp2));
                %mTh2 = full(median(gCCn((1:No)~=o, o)));
                %rAgrSum(j,k,o,l) = ( sum( (gCCn((1:No)~=o, o)>mTh2) & (tmp>mTh1)'  ) ...
                %    +  sum( (gCCn((1:No)~=o, o)<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
                rAgrSum(j,k,o,l) = ( sum( (tmp2>mTh2) & (tmp>mTh1)'  ) ...
                    +  sum( (tmp2<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
            end
        end
    end
    save(['pop4thOrderResp_iter_',num2str(popRespIter),'.mat'],'dSteps','kSteps','N4','r','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','Nsingle');
    
end

display({   '* * * *';'   * * * *';' calculation complete ';'   * * * *';'   * * * *'})

save('pop4thOrderResp.mat','dSteps','kSteps','N4','r','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','Nsingle');

% load('pop4thOrderResp.mat');
% rCCmean = mean(rCC,3);
% figure; plot( kSteps, rCCmean )
% hold on; plot(.1*ones(size(dSteps)), rCCmean(:,1), 'ko')
% %ylim([0 1])
% xlim([0 1])
% box off


rCCpMean = squeeze(mean(rCCpSum(:,1,:,:),3));
figure; plot( 1:N4, rCCpMean)
%hold on; plot(.1*ones(size(dSteps)), rCCpMean(:,1), 'ko')
%ylim([0 1])
xlim([0 N4])
xlabel('number of readout units')
box off


rCCmean = mean(rCC,3);
rCCsinglemean = mean(rCCSingle,3);

figure; semilogx( dSteps*10, rCCsinglemean )
hold all; semilogx( dSteps*10, rCCmean );
xlabel('number of piriform neurons')
legend('1 readout unit',[num2str(N4),' readout units'])
