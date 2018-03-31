

global Nx Ny Sx Sc Sp S4c NoR NoC classLevel1 classLevel2

% ------------
% ------------

% This just calculates parts for "untrained" 2 mouse fig 

% ------------
% ------------


tmpTime = clock;
rng(tmpTime(end)); %seed random number generator

initializeGlobals_Server;


NGA = Sx*Nx;
Cp = Sc*Nx; %glomerulus connections per piriform neuron
Cr = S4c*Ny; %piriform connections per 4th order neuron
trainedOdorNum = 1; %which odor to train (this is arbitrary)



% ODORS ---------------------------------------------------------

% first make correlation matrix:
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



% TRAINING & ABLATION ----------------------------------------------

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

tmp = rand(size(trainedWeights90p));
[~,loc] = sort(tmp,'descend');
UNTRAINEDWeights90p2 = trainedWeights90p(loc);

UNTRAINEDFourthOrder = UNTRAINEDWeights90p'*piriformRepBig;
UNTRAINEDFourthOrder_mouse2 = UNTRAINEDWeights90p2'*piriformRepBig_mouse2;

UNTRAINEDFourthOrder90pClass1 = UNTRAINEDFourthOrder(1:10*NoC/2);
UNTRAINEDFourthOrder90pClass2 = UNTRAINEDFourthOrder(10*NoC/2+1:10*NoC);
UNTRAINEDFourthOrder90pRand = UNTRAINEDFourthOrder(10*NoC+(1:10*NoC));

UNTRAINEDFourthOrder90pClass1_mouse2 = UNTRAINEDFourthOrder_mouse2(1:10*NoC/2);
UNTRAINEDFourthOrder90pClass2_mouse2 = UNTRAINEDFourthOrder_mouse2(10*NoC/2+1:10*NoC);
UNTRAINEDFourthOrder90pRand_mouse2 = UNTRAINEDFourthOrder_mouse2(10*NoC+(1:10*NoC));






%% SAVE EVERYTHING --------------------------------------------------------

save('partsForSumFig2.mat',...
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
    








