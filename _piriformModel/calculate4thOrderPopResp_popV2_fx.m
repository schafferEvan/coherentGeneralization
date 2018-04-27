
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
oCorr = 0;

dSteps = 1e6; %round(logspace(0,7,20)); %round(logspace(0,7,30));
pSteps = .5:.05:1; %.5:.05:.95; %[.7,.85,.95]; % fraction of required outputs needed to agree

%dSteps = dSteps(1:end-1);
%dSteps = dSteps(1:19);
%kSteps = 1;%.1:.2:.9;  %no longer in use
N4 = 1000; %50;
nperms = 2; %10;

glomRepLarge = makeOdors('correlated',oCorr,No); %makeOdors('correlated',0.7,No);

rAgr = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
rAgr90 = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
rAgrUn = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
rAgr90Un = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);

T50base = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
T90base = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
T50baseUn = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
T90baseUn = zeros(nperms, length(dSteps), N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
agr90base = .9^2+.1^2;


for j=1:length(dSteps)
    display(' ---------- ')
    display(j/length(dSteps))
    display(' ---------- ')
    
    Ny = dSteps(j);
    z = zeros(N4,No,No-1);
    zUn = zeros(N4,No,No-1);
    notOkFlag = zeros(N4,No);
    notOkFlagUn = zeros(N4,No);
    
    
    for k=1:N4
        display([j/length(dSteps), k/N4])
        
        % every z sees a different piriform
        piriformRepLarge = makePiriform(glomRepLarge); 
        zUnOrder = randperm(Ny); 
        
        % To minimize RAM pressure, compute each z and delete piriform
        for o=1:No
            
            % train readout to one odor
            vTot = -1/(No-1) * ones( N4, No );
            vTot(:,o) = 1;
            W = piriformRepLarge(:,o)';
            W = W-mean(mean(W));
            z(k,o,:) = full( W*piriformRepLarge(:,(1:No)~=o) ) / full( W*piriformRepLarge(:,o) );
            WUn = W(zUnOrder);
            zUn(k,o,:) = full( WUn*piriformRepLarge(:,(1:No)~=o) ) / full( WUn*piriformRepLarge(:,o) );
            
            if sum(sum(W))==0
                notOkFlag(k,o) = 1;
            end
            if sum(sum(WUn))==0
                notOkFlagUn(k,o) = 1;
            end
            
        end
        clear piriformRepLarge
    end
    
    T50base1 = repmat(median(z, 3),1,1,No-1);
    T90base1 = repmat(quantile(z, 0.9, 3),1,1,No-1);
    
    T50base1Un = repmat(median(zUn, 3),1,1,No-1);
    T90base1Un = repmat(quantile(zUn, 0.9, 3),1,1,No-1);
    
    % permutations of cumsum over N4
    for l=1:nperms
        pvec = randperm(N4);
        
        % for normalization, compute how many readouts have a valid weight
        % vector (nonzero W)
        ok2 = cumsum(~sum(isnan(z),3));
        ok2Un = cumsum(~sum(isnan(zUn),3));
        
        % fraction of readouts that respond/don't respond
        rAgrMat501 = cumsum(z(pvec,:,:)>T50base1,1).*repmat(1./ok2,1,1,No-1); 
        rAgrMat901 = cumsum(z(pvec,:,:)>T90base1,1).*repmat(1./ok2,1,1,No-1);
        %rAgrMat501 = cumsum(z(pvec,:,:)>T50base1,1).*repmat(1./(1:N4)',1,No,No-1); 
        %rAgrMat901 = cumsum(z(pvec,:,:)>T90base1,1).*repmat(1./(1:N4)',1,No,No-1);
        rAgrMat501Un = cumsum(zUn(pvec,:,:)>T50base1Un,1).*repmat(1./ok2Un,1,1,No-1); 
        rAgrMat901Un = cumsum(zUn(pvec,:,:)>T90base1Un,1).*repmat(1./ok2Un,1,1,No-1);
        
        
        for k=1:N4
            for o=1:No
                
                if notOkFlag(k,o)
                    rAgr(l,j,k,:,o) = nan;
                    rAgr90(l,j,k,:,o) = nan;
                    T50base(l,j,k,:,o) = nan;
                    T90base(l,j,k,:,o) = nan;
                    %continue
                else
                    for p=1:length(pSteps)
                        %defined fraction pSteps(p) required to count as agreement
                        rAgr(l,j,k,p,o) =   sum( (rAgrMat501(k,o,:)>pSteps(p)) | (rAgrMat501(k,o,:)<(1-pSteps(p))) )/(No-1); %sum( (rAgrMat501(k,o,:)>pSteps(p)) | (rAgrMat501(k,o,:)<1-pSteps(p)) )/(No-1);
                        rAgr90(l,j,k,p,o) = sum( (rAgrMat901(k,o,:)>pSteps(p)) | (rAgrMat901(k,o,:)<(1-pSteps(p))) )/(No-1);
                    end
                end
                
                if notOkFlagUn(k,o)
                    rAgrUn(l,j,k,:,o) = nan;
                    rAgr90Un(l,j,k,:,o) = nan;
                    T50baseUn(l,j,k,:,o) = nan;
                    T90baseUn(l,j,k,:,o) = nan;
                    %continue
                else
                    for p=1:length(pSteps)
                        %defined fraction pSteps(p) required to count as agreement
                        rAgrUn(l,j,k,p,o) = sum( rAgrMat501Un(k,o,:)>pSteps(p) )/(No-1); %sum( (rAgrMat501(k,o,:)>pSteps(p)) | (rAgrMat501(k,o,:)<1-pSteps(p)) )/(No-1);
                        %rAgr90(l,j,k,p,o) = sum( (rAgrMat901(k,o,:)>pSteps(p)) | (rAgrMat901(k,o,:)<1-pSteps(p)) )/(No-1);
                    end
                end

                
            end
        end
    end
    
    save(['pop4thOrderResp_POP_',num2str(oCorr),'_iter_',num2str(popRespIter),'_fx.mat'],'agr90base','T50base','T90base','rAgr90','oCorr','dSteps','pSteps','N4','rAgr');%,'T50baseUn','T90baseUn','rAgr90Un','rAgrUn');
    %figure; plot(squeeze(sum(sum(rAgrMat901,3),2)))
    
    %     for k=1:N4
    %         for o=1:No
    %
    %             if notOkFlag(k,o)
    %                 rAgr(j,k,:,o) = nan;
    %                 rAgr90(j,k,:,o) = nan;
    %                 T50base(j,k,:,o) = nan;
    %                 T90base(j,k,:,o) = nan;
    %                 continue
    %             end
    %
    %
    %
    %
    %             for p=1:length(pSteps)
    %                 %defined fraction pSteps(p) required to count as agreement
    %                 rAgr(j,k,p,o) = sum( (rAgrMat501(k,o,:)>pSteps(p)) | (rAgrMat501(k,o,:)<1-pSteps(p)) )/(No-1);
    %                 rAgr90(j,k,p,o) = sum( (rAgrMat901(k,o,:)>pSteps(p)) | (rAgrMat901(k,o,:)<1-pSteps(p)) )/(No-1);
    %             end
    %         end
    %     end
    %
    %     save(['pop4thOrderResp_POP_',num2str(oCorr),'_iter_',num2str(popRespIter),'.mat'],'agr90base','T50base','T90base','rAgr90','oCorr','dSteps','N4','rAgr');
    %
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
