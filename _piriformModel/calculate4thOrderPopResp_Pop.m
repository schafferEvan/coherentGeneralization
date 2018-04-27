
% updated WIdx for better memory management on 1/30/17

%% ...
global Ny

% for iterIter = 1:4; batchJob = batch(@calculate4thOrderPopResp_Pop,0);
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
No = 200; %800;%1000; %200;
oCorr = 0.7;

pSteps = [.5,.7,.9]; % fraction of required outputs needed to agree

N4 = 1000;
nInMax = 1000;

glomRepLarge = makeOdors('correlated',oCorr,No); %makeOdors('correlated',0.7,No);
gCC = corrcoef(glomRepLarge);

rAgr = zeros(N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);
rAgr90 = zeros(N4, length(pSteps), No); %zeros(Ny/dStep,dimIters);

agr90base = .9^2+.1^2;




gCCn = gCC;
for o=1:No
    gCCn((1:No)~=o, o) = gCCn((1:No)~=o, o) - mean(gCCn((1:No)~=o, o));
    gCCn((1:No)~=o, o) = gCCn((1:No)~=o, o) / norm(gCCn((1:No)~=o, o));
end

Ny = 1e5;
piriformRepLarge = makePiriform(glomRepLarge);
piriformRepLarge2 = makePiriform(glomRepLarge);




% restrict connectivity of readout
nIn = min(nInMax, size(piriformRepLarge,1));
%nIn = round(0.1* size(piriformRepLarge,1));


% Randomly select sparse connectivity for each 4th order neuron
WIdx = zeros(N4,nIn);
WIdx2 = zeros(N4,nIn);
for l=1:N4
    yOrder = randperm(Ny);
    WIdx(l, : ) = yOrder( 1:nIn );
    yOrder = randperm(Ny);
    WIdx2(l, : ) = yOrder( 1:nIn );
end


for o=1:No
    
    if ~mod(o/No,0.1)
        display(o/No)
    end
    
    % train all 4th order neurons to same odor
    vTot = -1/(No-1) * ones( N4, No );
    vTot(:,o) = 1;
    
    z = zeros(N4,No-1);
    z2 = zeros(N4,No-1);
    for jj=1:N4
        W = zeros(1,Ny);
        W2 = zeros(1,Ny);
        W(WIdx(jj,:)) = piriformRepLarge(WIdx(jj,:),o)';
        W2(WIdx2(jj,:)) = piriformRepLarge2(WIdx2(jj,:),o)';
        W = W-mean(mean(W));
        W2 = W2-mean(mean(W2));
        
        z(jj,:) = full( W*piriformRepLarge(:,(1:No)~=o) );
        z2(jj,:) = full( W2*piriformRepLarge2(:,(1:No)~=o) );
        %rvAll = full( W*piriformRepLarge ) .* vTot;
    end
    
    
    
        
    T50base1 = repmat(median(z, 2),1,No-1);
    T50base2 = repmat(median(z2, 2),1,No-1);
    T90base1 = repmat(quantile(z, 0.9, 2),1,No-1);
    T90base2 = repmat(quantile(z2, 0.9, 2),1,No-1);
    
    
    rAgrMat501 = cumsum(z>T50base1).*repmat(1./(1:N4)',1,No-1); % fraction of readouts that respond to each stimulus, based on chosen coding level
    rAgrMat502 = cumsum(z2>T50base2).*repmat(1./(1:N4)',1,No-1);
    rAgrMat901 = cumsum(z>T90base1).*repmat(1./(1:N4)',1,No-1);
    rAgrMat902 = cumsum(z2>T90base2).*repmat(1./(1:N4)',1,No-1);

    
    for p=1:length(pSteps)
        
        %% version 3
        % defined fraction pSteps(p) required to count as agreement
        %         rAgr(:,p,o) = ( sum( (rAgrMat502>pSteps(p)) & (rAgrMat501>pSteps(p)) , 2) ...
        %             +  sum( (rAgrMat502<=1-pSteps(p)) & (rAgrMat501<=1-pSteps(p)) , 2) )/(No-1);
        %         rAgr90(:,p,o) = ( sum( (rAgrMat902>pSteps(p)) & (rAgrMat901>pSteps(p)) , 2) ...
        %             +  sum( (rAgrMat902<=1-pSteps(p)) & (rAgrMat901<=1-pSteps(p)) , 2) )/(No-1);
        
        %% version 2
        % defined fraction pSteps(p) required to count as agreement
        rAgr(:,p,o) = ( sum( (rAgrMat502>pSteps(p)) & (rAgrMat501>pSteps(p)) , 2) ...
            +  sum( (rAgrMat502<=pSteps(p)) & (rAgrMat501<=pSteps(p)) , 2) )/(No-1);
        rAgr90(:,p,o) = ( sum( (rAgrMat902>pSteps(p)) & (rAgrMat901>pSteps(p)) , 2) ...
            +  sum( (rAgrMat902<=pSteps(p)) & (rAgrMat901<=pSteps(p)) , 2) )/(No-1);
        
        %% version 1
        %         rAgrMat = ( cumsum( (z2>T50base2) & (z>T50base1)  ) ...
        %             +  cumsum( (z2<T50base2) & (z<T50base1)  ) )...
        %             .*repmat(1./(1:N4)',1,No-1)/(No-1);
        %         rAgr90Mat = ( cumsum( (z2>T90base2) & (z>T90base1)  ) ...
        %             +  cumsum( (z2<T90base2) & (z<T90base1)  ) )...
        %             .*repmat(1./(1:N4)',1,No-1)/(No-1);
        %
        %         rAgr(:,k,o) = sum(rAgrMat,2);
        %         rAgr90(:,k,o) = sum(rAgr90Mat,2);

    end
    
end

save(['_pop4thOrderResp_0.7/_v2_forPop/pop4thOrderResp_POP_',num2str(oCorr),'_iter_',num2str(popRespIter),'.mat'],'agr90base','rAgr','rAgr90','oCorr','N4','pSteps');
display({   '* * * *';'   * * * *';' calculation complete ';'   * * * *';'   * * * *'})

figure; semilogx(1:N4,mean(rAgr,3))

%save('pop4thOrderResp.mat','dSteps','kSteps','N4','rAcc','rCC','gCC','rCCpSum','rCCSingle','rAgr','rAgrSingle','rAgrSum','RExSingle','REx','Nsingle','snr');
try
    cancel(batchMem)
catch
end
