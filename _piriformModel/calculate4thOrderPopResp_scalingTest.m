
%% temp plots
%figure; semilogx(dSteps,mean(squeeze(rAcc),2))
%hold all; semilogx(dSteps,mean(squeeze(rCC),2))
%hold all; semilogx(dSteps,-1+2*mean(squeeze(rAgr),2))
%hold all; semilogx(dSteps, mean(squeeze(snr),2)*( max(mean(squeeze(rCC),2))/max(mean(squeeze(snr),2))))


%% ...
global Ny

% for iterIter = 1:4; batchJob = batch(@calculate4thOrderPopResp_v1,0);
% wait(batchJob); end % doesn't work.  Needs input
%popRespIter = 0;


tmpTime = clock;
rng(1e5*tmpTime(end)); %seed random number generator
popRespIter = round(1e6*rand(1));
%popRespIter = 15; %popRespIter + 1;

if isunix && ~ismac % this actually would run on mac, but not wanted
    batchMem = batch(@logUnixMem,0,{popRespIter});
end

initializeGlobals_Server;
No = 100; %800;%1000; %200;
oCorr = 0;%0.7;

dSteps = round(logspace(1,4,30));%round(logspace(3,7,29)); %round(logspace(3,6.69,28)); %round(logspace(3,6,28)); %round(logspace(2,4,10)); %

%dSteps = dSteps(1:end-1);
%dSteps = dSteps(1:19);
kSteps = 1;%.1:.2:.9;  %no longer in use
N4 = 1; %200;

glomRepLarge = makeOdors('correlated',oCorr,No); %makeOdors('correlated',0.7,No);
%gCC = corrcoef(glomRepLarge);
rCC = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rCCpSum = zeros(length(dSteps), length(kSteps), No, N4); %zeros(Ny/dStep,dimIters);
rCCSingle = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rAcc = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);

rAgr = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rAgr90 = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);
rAgrSum = zeros(length(dSteps), length(kSteps), No, N4); %zeros(Ny/dStep,dimIters);
rAgrSingle = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);

T50base = zeros(length(dSteps), length(kSteps), No, 2); %zeros(Ny/dStep,dimIters);
T90base = zeros(length(dSteps), length(kSteps), No, 2); %zeros(Ny/dStep,dimIters);
agr90base = .9^2+.1^2;

snr = zeros(length(dSteps), length(kSteps), No); %zeros(Ny/dStep,dimIters);

yMu = zeros(length(dSteps), length(kSteps));
ySSq = zeros(length(dSteps), length(kSteps));

REx = zeros(length(dSteps), No-1);
REx2 = zeros(length(dSteps), No-1);
RExSingle = zeros(length(dSteps), No-1);
RExSingle2 = zeros(length(dSteps), No-1);
RExU = zeros(length(dSteps), No-1);
REx2U = zeros(length(dSteps), No-1);


for j=1:length(dSteps)
    display(j/length(dSteps))
    
    Ny = dSteps(j);
    piriformRepLarge = makePiriform(glomRepLarge);
    piriformRepLarge2 = makePiriform(glomRepLarge);
    
    for k=1:length(kSteps)
        
        % restrict connectivity of readout
        %nIn = min(1000, size(piriformRepLarge,1));
        %         nIn = Ny; %round(0.1* size(piriformRepLarge,1));
        
        %
        %         % Randomly select sparse connectivity for each 4th order neuron
        %         WIdx = zeros(N4,Ny);
        %         for l=1:N4
        %             yOrder = randperm(Ny);
        %             %WIdx(l, yOrder( 1:round(Ny*kSteps(k)) ) ) = 1;
        %             WIdx(l, yOrder( 1:nIn ) ) = 1;
        %         end
        %         WIdx2 = zeros(N4,Ny);
        %         for l=1:N4
        %             yOrder = randperm(Ny);
        %             %WIdx(l, yOrder( 1:round(Ny*kSteps(k)) ) ) = 1;
        %             WIdx2(l, yOrder( 1:nIn ) ) = 1;
        %         end
        
        yMu(j,k) = mean(reshape(piriformRepLarge,Ny*No,1));
        ySSq(j,k) = mean(var(piriformRepLarge,0,2));
            
        for o=1:No
            
            % train readout to one odor
            vTot = -1/(No-1) * ones( N4, No );
            vTot(:,o) = 1;
            W = piriformRepLarge(:,o)';
            W2 = piriformRepLarge2(:,o)';
            
            
            % new version post 12/20/16
            W = W-mean(mean(W));
            W2 = W2-mean(mean(W2));
            
            if sum(sum(W))==0 || sum(sum(W2))==0
                rAgr(j,k,o) = nan;
                rAgr90(j,k,o) = nan;
                rAcc(j,k,o) = nan;
                rCC(j,k,o) = nan;
                snr(j,k,o) = nan;
                %yMu(j,k,o) = nan;
                %ySSq(j,k,o) = nan;
                T50base(j,k,o) = nan;
                T90base(j,k,o) = nan;
                continue
            end
            
            % original version pre 12/20/16
            %W(W~=0) = W(W~=0) - mean(mean(W(W~=0)));
            %W2(W2~=0) = W2(W2~=0) - mean(mean(W2(W2~=0)));
            
            
            z = full( W*piriformRepLarge(:,(1:No)~=o) );
            z2 = full( W2*piriformRepLarge2(:,(1:No)~=o) );
            %rvAll = full( W*piriformRepLarge ) .* vTot;
            
            % Z is the "5th order" response, the sum of 4th order
            Z = sum(z,1);
            Z2 = sum(z2,1);
            
            
            
            tmp = Z-mean(Z);
            tmp = tmp/norm(tmp);
            tmp2 = Z2-mean(Z2);
            tmp2 = tmp2/norm(tmp2);
            %rCC(j,k,o) = tmp*gCCn((1:No)~=o, o); %tmp(2);
            rCC(j,k,o) = tmp*tmp2';
            
            REx(j,:) = tmp;
            REx2(j,:) = tmp2;
            mTh1 = full(median(tmp));
            mTh2 = full(median(tmp2));
            %mTh2 = full(median(gCCn((1:No)~=o, o)));
            %rAgr(j,k,o) = ( sum( (gCCn((1:No)~=o, o)>mTh2) & (tmp>mTh1)'  ) ...
            %    +  sum( (gCCn((1:No)~=o, o)<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
            T50base(j,k,o,1) = mTh1;
            T50base(j,k,o,2) = mTh2;
            T90base(j,k,o,1) = full(quantile(tmp,0.9));
            T90base(j,k,o,2) = full(quantile(tmp2,0.9));
            %rAgr(j,k,o) = ( sum( (tmp2>mTh2)' & (tmp>mTh1)'  ) ...
            %    +  sum( (tmp2<mTh2)' & (tmp<mTh1)'  ) )/length(tmp);
            %rAgr90(j,k,o) = ( sum( (tmp2>T90base(j,k,o,2))' & (tmp>T90base(j,k,o,1))'  ) ...
            %    +  sum( (tmp2<T90base(j,k,o,2))' & (tmp<T90base(j,k,o,1))'  ) )/length(tmp);
            valid50 = (tmp2~=mTh2) & (tmp~=mTh1);
            rAgr(j,k,o) = ( sum( (tmp2(valid50)>mTh2)' & (tmp(valid50)>mTh1)'  ) ...
                +  sum( (tmp2(valid50)<=mTh2)' & (tmp(valid50)<=mTh1)'  ) )/sum(valid50);
            valid90 = (tmp2~=T90base(j,k,o,2)) & (tmp~=T90base(j,k,o,1));
            rAgr90(j,k,o) = ( sum( (tmp2(valid90)>T90base(j,k,o,2))' & (tmp(valid90)>T90base(j,k,o,1))'  ) ...
                +  sum( (tmp2(valid90)<=T90base(j,k,o,2))' & (tmp(valid90)<=T90base(j,k,o,1))'  ) )/sum(valid90);
            
            %rAcc(j,k,o) = Z * vTot(1,(1:No)~=o)';
            rAcc(j,k,o) = sum( Z .* vTot(1,(1:No)~=o) < 0 ,2)/(No-1);
            
            
            %** intentionally backwards for "untrained"
            zU = full( W2*piriformRepLarge(:,(1:No)~=o) );
            z2U = full( W*piriformRepLarge2(:,(1:No)~=o) );            
            % Z is the "5th order" response, the sum of 4th order
            ZU = sum(zU,1);
            Z2U = sum(z2U,1);
            tmpU = ZU-mean(ZU);
            tmpU = tmpU/norm(tmpU);
            tmp2U = Z2U-mean(Z2U);
            tmp2U = tmp2U/norm(tmp2U);            
            RExU(j,:) = tmpU;
            REx2U(j,:) = tmp2U;

            
            Nsingle = 10;
            tmpV = zeros(1,floor(N4/Nsingle));
            tmpV2 = zeros(1,floor(N4/Nsingle));
            for l=1:length(tmpV)
                tmp = sum(z( 1+(l-1)*Nsingle:l*Nsingle, :));
                tmp2 = sum(z2( 1+(l-1)*Nsingle:l*Nsingle, :));
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
                tmpV2 = ( sum( (tmp2>mTh2)' & (tmp>mTh1)'  ) ...
                    +  sum( (tmp2<mTh2)' & (tmp<mTh1)'  ) )/length(tmp);
            end
            %rCCSingle(j,k,o) = mean(tmpV);
            rAgrSingle(j,k,o) = mean(tmpV2);
            RExSingle(j,:) = tmp;
            
            %r(j,k,:,:) = z'; %this is wrong
            
            Zc = cumsum(z);
            Zc2 = cumsum(z2);
            
            % compute mahalanobis distance between trained odor & others
            %trainedResp = full( W*piriformRepLarge(:,o) );
            %trainedResp2 = full( W2*piriformRepLarge2(:,o) );
            %T = [sum(trainedResp),sum(trainedResp2)];
            %snr(j,k,o) = sqrt( mahal(T, [Zc(end,:)',Zc2(end,:)']) );
            %snr(j,k,o) = sum(sum(rvAll))^2 / var(reshape(rvAll,N4*No,1));
            snr(j,k,o) = mean(mean(z))^2 / var(reshape(z,N4*(No-1),1));
            
            
            for l=1:N4
                tmp = Zc(l,:)-mean(Zc(l,:));
                tmp = tmp/norm(tmp);
                tmp2 = Zc2(l,:)-mean(Zc2(l,:));
                tmp2 = tmp2/norm(tmp2);
                %tmp = corrcoef( Zc(l,:), gCC((1:No)~=o, o) );
                %rCCpSum(j,k,o,l) = tmp*gCCn((1:No)~=o, o); %tmp(2);
                
                mTh1 = full(median(tmp));
                mTh2 = full(median(tmp2));
                %mTh2 = full(median(gCCn((1:No)~=o, o)));
                %rAgrSum(j,k,o,l) = ( sum( (gCCn((1:No)~=o, o)>mTh2) & (tmp>mTh1)'  ) ...
                %    +  sum( (gCCn((1:No)~=o, o)<mTh2) & (tmp<mTh1)'  ) )/length(tmp);
                rAgrSum(j,k,o,l) = ( sum( (tmp2>mTh2)' & (tmp>mTh1)'  ) ...
                    +  sum( (tmp2<mTh2)' & (tmp<mTh1)'  ) )/length(tmp);
            end
        end
    end
    clear piriformRepLarge piriformRepLarge2
    save(['pop4thOrderResp_scaleTest_iter_',num2str(popRespIter),'.mat'],'yMu','ySSq','agr90base','T50base','T90base','rAgr90','oCorr','dSteps','kSteps','N4','rAcc','rCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','REx2','RExU','REx2U','snr');
    
end

display({   '* * * *';'   * * * *';' calculation complete ';'   * * * *';'   * * * *'})

%save('pop4thOrderResp.mat','dSteps','kSteps','N4','rAcc','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','Nsingle','snr');
try
    cancel(batchMem)
catch
end

% load('pop4thOrderResp.mat');
% rCCmean = mean(rCC,3);
% figure; plot( kSteps, rCCmean )
% hold on; plot(.1*ones(size(dSteps)), rCCmean(:,1), 'ko')
% %ylim([0 1])
% xlim([0 1])
% box off
% 
% 
% rCCpMean = squeeze(mean(rCCpSum(:,1,:,:),3));
% figure; plot( 1:N4, rCCpMean)
% %hold on; plot(.1*ones(size(dSteps)), rCCpMean(:,1), 'ko')
% %ylim([0 1])
% xlim([0 N4])
% xlabel('number of readout units')
% box off
% 
% 
% rCCmean = mean(rCC,3);
% rCCsinglemean = mean(rCCSingle,3);
% 
% figure; semilogx( dSteps*10, rCCsinglemean )
% hold all; semilogx( dSteps*10, rCCmean );
% xlabel('number of piriform neurons')
% legend('1 readout unit',[num2str(N4),' readout units'])
