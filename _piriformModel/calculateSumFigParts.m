

global Nx Ny Sx S4c Sc Sp NoR NoC classLevel1 classLevel2


tmpTime = clock;
rng(1e5*tmpTime(end)); %seed random number generator

initializeGlobals_Server;

% updated to use only "correlated" odors instead of mix & class on 12/29/15
% split into plotting part and data generating part (calculateSumFigParts.m) on 1/5/16


NGA = Sx*Nx;
Cp = Sc*Nx; %glomerulus connections per piriform neuron
Cr = S4c*Ny; %piriform connections per 4th order neuron
trainedOdorNum = 1; %which odor to train (this is arbitrary)



% ODORS ---------------------------------------------------------
NoTot = NoR+NoC;
Sclass = classLevel2*(ones(NoC)-eye(NoC))+eye(NoC);
Sclass(1:NoC/2,1:NoC/2) = classLevel1*(ones(NoC/2)-eye(NoC/2))+eye(NoC/2);
S = eye(NoTot);
S(1:NoC,1:NoC) = Sclass;


[glomRepTemp,odorSimilarity] = makeOdors('correlated',S,NoTot);
glomRepTempMix = nan(Nx,NoM);
odorSimilarity2 = nan(1,NoM);
ds = 50;
idx = 1+(0:ds:NoM);

for j=1:length(idx)-1
    [glomRepTempMix(:,idx(j):idx(j+1)-1), odorSimilarity2(:,idx(j):idx(j+1)-1)] = makeOdors('correlated',j/(length(idx)-1),ds); %note: these mini-classes can't be compared to each other (not created simultaneously)
end

similarityMeasureFull = [odorSimilarity, -odorSimilarity2];
similarityMat = repmat(similarityMeasureFull',1,length(similarityMeasureFull)) == repmat(similarityMeasureFull,length(similarityMeasureFull),1);
glomerulusRep = [glomRepTemp, glomRepTempMix];

% PIRIFORM -------------------------------------------------------

[piriformRep,~,~,J,maxy] = makePiriform(glomerulusRep);  % mouse 1
[~,~,~,J2] = makePiriform(glomerulusRep,[],[],maxy);  % mouse 1

% CORRELATIONS & DISTANCES
glomCorr = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
pirCorr = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
pirCorr10xSparse = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
pirCorr100xSparse = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
pirCorr1000xSparse = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
glomOverlap = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
pirOverlap = NaN*ones((NoR+NoC+NoM)); %NaN*ones((NoR+NoC+NoM)/2,(NoR+NoC+NoM));
glomDist = NaN*ones((NoR+NoC+NoM));
pirDist = NaN*ones((NoR+NoC+NoM));

corrIdx1 = NaN*ones((NoR+NoC+NoM));
corrIdx2 = NaN*ones((NoR+NoC+NoM));

glomRepNorm = glomerulusRep;
pirRepNorm = piriformRep;
[~,loc] = sort(rand(Ny,1),'descend');
pirRepNorm10xSparse = piriformRep(loc(1:round(.1*Ny)),:);
pirRepNorm100xSparse = piriformRep(loc(1:round(.01*Ny)),:);
pirRepNorm1000xSparse = piriformRep(loc(1:round(.001*Ny)),:);

for k=1:(NoR+NoC+NoM)
    glomRepNorm(:,k) = glomRepNorm(:,k) - mean(glomRepNorm(:,k));
    glomRepNorm(:,k) = glomRepNorm(:,k)/norm(glomRepNorm(:,k));
    pirRepNorm(:,k) = pirRepNorm(:,k) - mean(pirRepNorm(:,k));
    pirRepNorm(:,k) = pirRepNorm(:,k)/norm(pirRepNorm(:,k));
    pirRepNorm10xSparse(:,k) = pirRepNorm10xSparse(:,k) - mean(pirRepNorm10xSparse(:,k));
    pirRepNorm10xSparse(:,k) = pirRepNorm10xSparse(:,k)/norm(pirRepNorm10xSparse(:,k));
    pirRepNorm100xSparse(:,k) = pirRepNorm100xSparse(:,k) - mean(pirRepNorm100xSparse(:,k));
    pirRepNorm100xSparse(:,k) = pirRepNorm100xSparse(:,k)/norm(pirRepNorm100xSparse(:,k));
    pirRepNorm1000xSparse(:,k) = pirRepNorm1000xSparse(:,k) - mean(pirRepNorm1000xSparse(:,k));
    pirRepNorm1000xSparse(:,k) = pirRepNorm1000xSparse(:,k)/norm(pirRepNorm1000xSparse(:,k));
end


for i=1:(NoR+NoC+NoM)
    for ii=1:i-1
        
        if ~similarityMat(i,ii)
            continue; % only compare odors of the same type
        end
        
        glomCorr(i,ii) = glomRepNorm(:,i)'*glomRepNorm(:,ii); %glomtemp(2);
        pirCorr(i,ii) = pirRepNorm(:,i)'*pirRepNorm(:,ii); %pirtemp(2);
        pirCorr10xSparse(i,ii) = pirRepNorm10xSparse(:,i)'*pirRepNorm10xSparse(:,ii); %pirtemp(2);
        pirCorr100xSparse(i,ii) = pirRepNorm100xSparse(:,i)'*pirRepNorm100xSparse(:,ii); %pirtemp(2);
        pirCorr1000xSparse(i,ii) = pirRepNorm1000xSparse(:,i)'*pirRepNorm1000xSparse(:,ii); %pirtemp(2);
        corrIdx1(i,ii) = i;
        corrIdx2(i,ii) = ii;
         
        glomOverlap(i,ii) = 1/(Nx)*(glomerulusRep(:,i)>0)'*(glomerulusRep(:,ii)>0);
        pirOverlap(i,ii) = 1/(Ny)*(piriformRep(:,i)>0)'*(piriformRep(:,ii)>0);
        
        glomDist(i,ii) = (1/Nx^2)*sum( ( glomerulusRep(:,i) - glomerulusRep(:,ii) ).^2 );
        pirDist(i,ii) = (1/Ny^2)*sum( ( piriformRep(:,i) - piriformRep(:,ii) ).^2 );
        
    end
end

glomCorr = glomCorr(~isnan(glomCorr));
pirCorr = pirCorr(~isnan(pirCorr));
pirCorr10xSparse = pirCorr10xSparse(~isnan(pirCorr10xSparse));
pirCorr100xSparse = pirCorr100xSparse(~isnan(pirCorr100xSparse));
pirCorr1000xSparse = pirCorr1000xSparse(~isnan(pirCorr1000xSparse));

glomDist = glomDist(~isnan(glomDist));
pirDist = pirDist(~isnan(pirDist));

corrIdx1 = corrIdx1(~isnan(corrIdx1));
corrIdx2 = corrIdx2(~isnan(corrIdx2));
corrIsClass1Idx = (corrIdx1<=NoC/2) & (corrIdx2<=NoC/2);
corrIsClass2Idx = (corrIdx1<=NoC) & (corrIdx2<=NoC) & (corrIdx1>NoC/2) & (corrIdx2>NoC/2);
corrIsRandIdx = (corrIdx1>NoC) & (corrIdx2>NoC) & (corrIdx1<=(NoR+NoC)) & (corrIdx2<=(NoR+NoC));
corrIsMixtureIdx = (corrIdx1>(NoR+NoC)) & (corrIdx2>(NoR+NoC));

% first make correlation matrix:
NoTot = 10*(NoR+NoC);
Sclass = classLevel2*(ones(10*NoC)-eye(10*NoC))+eye(10*NoC);
Sclass(1:10*NoC/2,1:10*NoC/2) = classLevel1*(ones(10*NoC/2)-eye(10*NoC/2))+eye(10*NoC/2);
S = eye(NoTot);
S(1:10*NoC,1:10*NoC) = Sclass;


glomRepTempExp = makeOdors('correlated',S,NoTot);

piriformRepBig = makePiriform(glomRepTempExp,Sp,J,maxy);  % mouse 1

piriformRepBig_mouse2 = makePiriform(glomRepTempExp,Sp,J2,maxy); % mouse 2

s=whos('piriformRepBig_mouse2');
memBytes = s.bytes;
display(memBytes)

trainedOdor = piriformRepBig(:,trainedOdorNum);    % mouse 1
trainedOdor2 = piriformRepBig_mouse2(:,trainedOdorNum);  % mouse 2

trainedWeights = trainedOdor; %(trainedOdor~=0);
trainedWeights2 = trainedOdor2; %(trainedOdor2~=0);

nonzeroInputs = find(trainedWeights);
nonzeroInputs2 = find(trainedWeights2);

numUsedWeights = sum(trainedOdor~=0); %Note: this can be Sp*Ny if # cells/odor is exact rather than variable
numUsedWeights2 = sum(trainedOdor2~=0);

[~,ablationloc] = sort(randn(numUsedWeights,1),'descend');
[~,ablationloc2] = sort(randn(numUsedWeights2,1),'descend');

trainedWeights90p = trainedWeights;
trainedWeights90p(nonzeroInputs(ablationloc(1:round(numUsedWeights*.9)))) = 0;
trainedWeights90p2 = trainedWeights2;
trainedWeights90p2(nonzeroInputs2(ablationloc2(1:round(numUsedWeights2*.9)))) = 0;

trainedWeights90p = trainedWeights90p - mean(trainedWeights90p);
trainedWeights90p2 = trainedWeights90p2 - mean(trainedWeights90p2);


trainedFourthOrder = trainedWeights90p'*piriformRepBig;
trainedFourthOrder_mouse2 = trainedWeights90p2'*piriformRepBig_mouse2;

trainedFourthOrder90pClass1 = trainedFourthOrder(1:10*NoC/2);
trainedFourthOrder90pClass2 = trainedFourthOrder(10*NoC/2+1:10*NoC);
trainedFourthOrder90pRand = trainedFourthOrder(10*NoC+(1:10*NoC));

trainedFourthOrder90pClass1_mouse2 = trainedFourthOrder_mouse2(1:10*NoC/2);
trainedFourthOrder90pClass2_mouse2 = trainedFourthOrder_mouse2(10*NoC/2+1:10*NoC);
trainedFourthOrder90pRand_mouse2 = trainedFourthOrder_mouse2(10*NoC+(1:10*NoC));

% UNTRAINED READOUT: this is just swapped & shuffled readout to achieve same stats
tmp = rand(size(trainedWeights90p));
[~,loc] = sort(tmp,'descend');
UNTRAINEDWeights90p = trainedWeights90p2(loc);

%tmp = rand(size(trainedWeights90p));
%[~,loc] = sort(tmp,'descend');
UNTRAINEDWeights90p2 = UNTRAINEDWeights90p; %trainedWeights90p(loc);

UNTRAINEDFourthOrder = UNTRAINEDWeights90p'*piriformRepBig;
UNTRAINEDFourthOrder_mouse2 = UNTRAINEDWeights90p2'*piriformRepBig_mouse2;

UNTRAINEDFourthOrder90pClass1 = UNTRAINEDFourthOrder(1:10*NoC/2);
UNTRAINEDFourthOrder90pClass2 = UNTRAINEDFourthOrder(10*NoC/2+1:10*NoC);
UNTRAINEDFourthOrder90pRand = UNTRAINEDFourthOrder(10*NoC+(1:10*NoC));

UNTRAINEDFourthOrder90pClass1_mouse2 = UNTRAINEDFourthOrder_mouse2(1:10*NoC/2);
UNTRAINEDFourthOrder90pClass2_mouse2 = UNTRAINEDFourthOrder_mouse2(10*NoC/2+1:10*NoC);
UNTRAINEDFourthOrder90pRand_mouse2 = UNTRAINEDFourthOrder_mouse2(10*NoC+(1:10*NoC));





%% SAVE EVERYTHING --------------------------------------------------------

save('partsForSumFig.mat',...
    'glomCorr',...
    'corrIsRandIdx',...
    'corrIsClass1Idx',...
    'corrIsClass2Idx',...
    'corrIsMixtureIdx',...
    'pirCorr',...
    'pirCorr10xSparse',...
    'pirCorr100xSparse',...
    'pirCorr1000xSparse',...
    'glomOverlap',...
    'pirOverlap',...
    'glomDist',...
    'pirDist',...
    'trainedFourthOrder90pClass1',...
    'trainedFourthOrder90pClass2',...
    'trainedFourthOrder90pRand',...
    'trainedFourthOrder90pClass1_mouse2',...
    'trainedFourthOrder90pClass2_mouse2',...
    'trainedFourthOrder90pRand_mouse2',...
    'trainedOdorNum',...
    'UNTRAINEDFourthOrder90pClass1',...
    'UNTRAINEDFourthOrder90pClass2',...
    'UNTRAINEDFourthOrder90pRand',...
    'UNTRAINEDFourthOrder90pClass1_mouse2',...
    'UNTRAINEDFourthOrder90pClass2_mouse2',...
    'UNTRAINEDFourthOrder90pRand_mouse2',...
    'memBytes');
    








