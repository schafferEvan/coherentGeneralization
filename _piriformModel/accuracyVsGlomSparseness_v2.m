
% for reviewer

%% ...
global Ny Sx thMean;

tmpTime = clock;
rng(1e5*tmpTime(end)); %seed random number generator
popRespIter = round(1e6*rand(1));

initializeGlobals_Server;
No = 100; %800;%1000; %200;
oCorr = 0;%0.7;

dSteps = 100000;
xsSteps = .02:.02:.3; %input sparseness
N4 = 1; %200;
iters = 1; %10;

th1 = thMean;
th2 = thMean;

rCC = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);
rCCpSum = zeros(length(dSteps), length(xsSteps), No, N4, iters); %zeros(Ny/dStep,dimIters);
rCCSingle = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);
rAcc = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);

rAgr = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);
rAgr90 = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);
rAgrSum = zeros(length(dSteps), length(xsSteps), No, N4, iters); %zeros(Ny/dStep,dimIters);
rAgrSingle = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);

T50base = zeros(length(dSteps), length(xsSteps), No, 2, iters); %zeros(Ny/dStep,dimIters);
T90base = zeros(length(dSteps), length(xsSteps), No, 2, iters); %zeros(Ny/dStep,dimIters);
agr90base = .9^2+.1^2;

snr = zeros(length(dSteps), length(xsSteps), No, iters); %zeros(Ny/dStep,dimIters);

ZCat = zeros(length(dSteps), length(xsSteps), No*(No-1));
Z2Cat = zeros(length(dSteps), length(xsSteps), No*(No-1));

REx = zeros(length(dSteps), No-1);
REx2 = zeros(length(dSteps), No-1);
RExSingle = zeros(length(dSteps), No-1);
RExSingle2 = zeros(length(dSteps), No-1);
RExU = zeros(length(dSteps), No-1);
REx2U = zeros(length(dSteps), No-1);

for i=1:iters
    for j=1:length(dSteps)
        for k=1:length(xsSteps)
            display([j/length(dSteps),k/length(xsSteps)])
            
            Sx = xsSteps(k);
            glomRepLarge = makeOdors('correlated',oCorr,No); %makeOdors('correlated',0.7,No);
            
            Ny = dSteps(j);
            %[piriformRepLarge,~,~,~,~,th1] = makePiriform(glomRepLarge,[],[],[],[],th1);
            %[piriformRepLarge2,~,~,~,~,th2] = makePiriform(glomRepLarge,[],[],[],[],th2);
            piriformRepLarge = makePiriform(glomRepLarge,[],[],[],[],th1);
            piriformRepLarge2 = makePiriform(glomRepLarge,[],[],[],[],th2); %enforces fixed threshold

            
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
                    ZCat(j,k,(1:No-1)+(No-1)*(o-1)) = nan;
                    Z2Cat(j,k,(1:No-1)+(No-1)*(o-1)) = nan;
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
                
                ZCat(j,k,(1:No-1)+(No-1)*(o-1)) = Z;
                Z2Cat(j,k,(1:No-1)+(No-1)*(o-1)) = Z2;
                
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
                T50base(j,k,o,1,i) = mTh1;
                T50base(j,k,o,2,i) = mTh2;
                T90base(j,k,o,1,i) = full(quantile(tmp,0.9));
                T90base(j,k,o,2,i) = full(quantile(tmp2,0.9));
                %rAgr(j,k,o) = ( sum( (tmp2>mTh2)' & (tmp>mTh1)'  ) ...
                %    +  sum( (tmp2<mTh2)' & (tmp<mTh1)'  ) )/length(tmp);
                %rAgr90(j,k,o) = ( sum( (tmp2>T90base(j,k,o,2))' & (tmp>T90base(j,k,o,1))'  ) ...
                %    +  sum( (tmp2<T90base(j,k,o,2))' & (tmp<T90base(j,k,o,1))'  ) )/length(tmp);
                valid50 = (tmp2~=mTh2) & (tmp~=mTh1);
                rAgr(j,k,o,i) = ( sum( (tmp2(valid50)>mTh2)' & (tmp(valid50)>mTh1)'  ) ...
                    +  sum( (tmp2(valid50)<=mTh2)' & (tmp(valid50)<=mTh1)'  ) )/sum(valid50);
                valid90 = (tmp2~=T90base(j,k,o,2)) & (tmp~=T90base(j,k,o,1));
                rAgr90(j,k,o,i) = ( sum( (tmp2(valid90)>T90base(j,k,o,2))' & (tmp(valid90)>T90base(j,k,o,1))'  ) ...
                    +  sum( (tmp2(valid90)<=T90base(j,k,o,2))' & (tmp(valid90)<=T90base(j,k,o,1))'  ) )/sum(valid90);
                
                %rAcc(j,k,o) = Z * vTot(1,(1:No)~=o)';
                rAcc(j,k,o,i) = sum( Z .* vTot(1,(1:No)~=o) < 0 ,2)/(No-1);
                
                
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
                rAgrSingle(j,k,o,i) = mean(tmpV2);
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
                snr(j,k,o,i) = sum(sum(z))^2 / var(reshape(z,N4*(No-1),1));
                
                
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
                    rAgrSum(j,k,o,l,i) = ( sum( (tmp2>mTh2)' & (tmp>mTh1)'  ) ...
                        +  sum( (tmp2<mTh2)' & (tmp<mTh1)'  ) )/length(tmp);
                end
            end
        end
        clear piriformRepLarge piriformRepLarge2
        save(['popRespVsSparsenss_V2_',num2str(oCorr),'_iter_',num2str(popRespIter),'.mat'],'ZCat','Z2Cat','agr90base','T50base','T90base','rAgr90','oCorr','dSteps','xsSteps','N4','rAcc','rCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','REx2','RExU','REx2U','snr');
        
    end
end
disp({   '* * * *';'   * * * *';' calculation complete ';'   * * * *';'   * * * *'})

makePlots = false;
if makePlots
    figure;
    set(gcf,'color','w')
    set(gcf,'Units','Centimeters')
    pos(3)=8; %11.4; %17.4; %8.5;%pos(3)*1.5;
    pos(4)=8; %14.0;%pos(4)*1.5;
    set(gcf,'Position',pos)
    plot( xsSteps, mean(mean(snr,4),3)/max(mean(mean(snr,4),3)), 'linewidth',2, 'color',[.8,.8,.8]); hold on
    plot( xsSteps, mean(mean(rCC,4),3)/max(mean(mean(rCC,4),3)), '--','linewidth',1, 'color',[.2,.2,.2]); hold on
    plot( xsSteps, mean(mean(rAcc,4),3)/max(mean(mean(rAcc,4),3)), 'linewidth',1, 'color',[.8 0 .4] );
    rAgrmeanNorm = mean(mean(rAgr,4),3);
    rAgrmeanNorm = 2*(rAgrmeanNorm - 0.5);
    plot( xsSteps, rAgrmeanNorm/max(rAgrmeanNorm), 'linewidth',1, 'color',[.6 .3 0] );
    plot( xsSteps, 1/(1-agr90base)*(-agr90base+mean(mean(rAgr90,4),3)/max(mean(mean(rAgr90,4),3))), '--', 'linewidth',1, 'color',[.9 .6 0] )
    box off
    xlabel('Input Sparseness','fontsize',11)
    ylabel('Fraction of Maximum','fontsize',11)
    [hl,legHndl] = legend('SNR','Correlation','Accuracy','Agreement (\theta=0.5)','Agreement (\theta=0.9)','Location','best'); %'Southeast');
    %[hl,legHndl] = legend('Correlation','Accuracy','Agreement (\theta=0.5)','Agreement (\theta=0.9)','Location','best'); %'Southeast');
    set(hl,'FontSize',10);
    legend boxoff
    hLines=findobj(legHndl,'type','line');
    set(hLines,'linewidth',1.5);
end
