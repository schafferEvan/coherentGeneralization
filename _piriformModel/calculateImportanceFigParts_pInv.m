
%function makeImportanceFig
global  Ny Sp NoR NoC

%rng(second(datetime('now'))); %seed random number generator
tmpTime = clock;
rng(tmpTime(end)); %seed random number generator

initializeGlobals_Server;
%initializeGlobals;

% changed to "correlated" odors on 1/5/16
% split into plotting part and data generating part (calculateImportanceFigParts.m) on 1/5/16
% reorganized to reduce RAM load on 1/6/16

numOdorsToPlot = 200;
maxy = 17.3607; %this is to avoid massive memory load to find maxy on large dataset

NoR6 = NoR*5; % NOTE: NoR6 and NoC6 must be the same
NoC6 = NoC*5;
NoM6 = 100;
normalizeResp = false;

% ODORS ---------------------------------------------------------

[glomRepTempRand,odorSimilarity] = makeOdors('correlated',0,NoR6);
similarityMeasureFull = -eps*(odorSimilarity<0);
[glomRepTempClass,odorSimilarity] = makeOdors('correlated',classLevel1,NoC6); 
similarityMeasureFull = [-2+eps*(odorSimilarity<0),similarityMeasureFull];

glomRepMouse2 = makeOdors('correlated',0,NoR6);

% PIRIFORM -------------------------------------------------------

[pirClass,~,~,J] = makePiriform(glomRepTempClass,Sp,[],maxy);  % mouse 1
pirRand = makePiriform(glomRepTempRand,Sp,J,maxy);  % mouse 1
pirClassRand = makePiriform(glomRepMouse2,Sp,J,maxy); % mouse 2


%%  CLASS DISCRIMINANT & IMPORTANCE (class example) -------------------------
mc = mean(pirClass,2);
mr = mean(pirRand,2);

importance = (mc-mr);% ./ mr;% ./ (mc+mr);

[~,importanceLoc] = sort(importance,'descend');
importanceBound1a = round(.03*Ny);
importanceBound1b = round(.04*Ny);
importanceBound2 = find(importance(importanceLoc)>0, 1, 'last'); %Ny;
importanceSparseBound = 0.1;%Sc; %.99;
[~,randLoc] = sort(rand(size(mc)),'descend');
sparseLocVec = zeros(size(mc));
sparseLocVec(randLoc(1:round(Ny*importanceSparseBound))) = 1;

A = full([pirClass,pirRand]);
%mtemp = mean(mean(A));
classWeightsFull = [ones(1,NoC6),-ones(1,NoR6)]*pinv(A);%-mtemp);
Apart = A(importanceLoc(importanceBound1a+1:importanceBound2),:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsNoTop = zeros(size(classWeightsFull));
classWeightsNoTop(importanceLoc(importanceBound1a+1:importanceBound2)) = tempW;
Apart = A(importanceLoc(1:importanceBound1a),:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsNoBottom = zeros(size(classWeightsFull));
classWeightsNoBottom(importanceLoc(1:importanceBound1a)) = tempW;
Apart = A(sparseLocVec>0,:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsRand = zeros(size(classWeightsFull));
classWeightsRand(sparseLocVec>0) = tempW;


noiseVec = .05*randn(1,NoC6+NoR6);
jitterMat = std(A,[],2)*noiseVec;
At = A+jitterMat;


classFourthOrderFull = full(classWeightsFull*At);
classFourthOrderNoTop = full(classWeightsNoTop*At);
classFourthOrderNoBottom = full(classWeightsNoBottom*At);
classFourthOrderRand = full(classWeightsRand*At);

if normalizeResp
    %normalization & shifting
    classFourthOrderFull = (classFourthOrderFull-min(classFourthOrderFull))/(max(classFourthOrderFull)-min(classFourthOrderFull));
    classFourthOrderNoTop = (classFourthOrderNoTop-min(classFourthOrderNoTop))/(max(classFourthOrderNoTop)-min(classFourthOrderNoTop));
    classFourthOrderNoBottom = (classFourthOrderNoBottom-min(classFourthOrderNoBottom))/(max(classFourthOrderNoBottom)-min(classFourthOrderNoBottom));
    classFourthOrderRand = (classFourthOrderRand-min(classFourthOrderRand))/(max(classFourthOrderRand)-min(classFourthOrderRand));
else
    %just shifting
    classFourthOrderFull = classFourthOrderFull-min(classFourthOrderFull);
    classFourthOrderNoTop = classFourthOrderNoTop-min(classFourthOrderNoTop);
    classFourthOrderNoBottom = classFourthOrderNoBottom-min(classFourthOrderNoBottom);
    classFourthOrderRand = classFourthOrderRand-min(classFourthOrderRand);
end

selectedNeuron = importanceLoc(importanceBound1a:importanceBound1b);
classResponsesPt1 = pirClass(selectedNeuron,1:numOdorsToPlot);
randResponsesPt1 = pirRand(selectedNeuron,1:numOdorsToPlot);

save('part1ForImportanceFig_pInv.mat','selectedNeuron','classResponsesPt1','randResponsesPt1');
clear pirClass
    
%% CLASS DISCRIMINANT & IMPORTANCE (RANDOM example) -------------------------
mcr = mean(pirClassRand,2);
importanceRand = (mcr-mr);% ./ mr;% ./ (mcr+mr);

[~,importanceLocRand] = sort(importanceRand,'descend');
importanceBound1aRand = round(.03*Ny);
importanceBound1bRand = round(.04*Ny);
importanceBound2Rand = find(importanceRand(importanceLocRand)>0, 1, 'last'); %Ny;
temp = rand(size(mcr))<importanceSparseBound;

A = full([pirClassRand,pirRand]);
mtemp = mean(mean(A));
classWeightsFullRand = [ones(1,NoR6),-ones(1,NoR6)]*pinv(A); %-mtemp);
Apart = A(importanceLocRand(importanceBound1aRand+1:importanceBound2Rand),:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsNoTopRand = zeros(size(classWeightsFullRand));
classWeightsNoTopRand(importanceLocRand(importanceBound1aRand+1:importanceBound2Rand)) = tempW;
Apart = A(importanceLocRand(1:importanceBound1aRand),:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsNoBottomRand = zeros(size(classWeightsFullRand));
classWeightsNoBottomRand(importanceLocRand(1:importanceBound1aRand)) = tempW;
Apart = A(sparseLocVec>0,:);
tempW = [ones(1,NoC6),-ones(1,NoR6)]*pinv(Apart);
classWeightsRandRand = zeros(size(classWeightsFullRand));
classWeightsRandRand(sparseLocVec>0) = tempW;


noiseVec = .05*randn(1,NoC6+NoR6);
jitterMat = std(A,[],2)*noiseVec;
At = A+jitterMat;
classFourthOrderFullRand = full(classWeightsFullRand*At);
classFourthOrderNoTopRand = full(classWeightsNoTopRand*At);
classFourthOrderNoBottomRand = full(classWeightsNoBottomRand*At);
classFourthOrderRandRand = full(classWeightsRandRand*At);

if normalizeResp
    %normalization & shifting
    classFourthOrderFullRand = (classFourthOrderFullRand-min(classFourthOrderFullRand))/(max(classFourthOrderFullRand)-min(classFourthOrderFullRand));
    classFourthOrderNoTopRand = (classFourthOrderNoTopRand-min(classFourthOrderNoTopRand))/(max(classFourthOrderNoTopRand)-min(classFourthOrderNoTopRand));
    classFourthOrderNoBottomRand = (classFourthOrderNoBottomRand-min(classFourthOrderNoBottomRand))/(max(classFourthOrderNoBottomRand)-min(classFourthOrderNoBottomRand));
    classFourthOrderRandRand = (classFourthOrderRandRand-min(classFourthOrderRandRand))/(max(classFourthOrderRandRand)-min(classFourthOrderRandRand));
else
    %just shifting
    classFourthOrderFullRand = classFourthOrderFullRand-min(classFourthOrderFullRand);
    classFourthOrderNoTopRand = classFourthOrderNoTopRand-min(classFourthOrderNoTopRand);
    classFourthOrderNoBottomRand = classFourthOrderNoBottomRand-min(classFourthOrderNoBottomRand);
    classFourthOrderRandRand = classFourthOrderRandRand-min(classFourthOrderRandRand);
end

selectedNeuronRand = importanceLocRand(importanceBound1aRand:importanceBound1bRand);
randResponsesPt2 = pirRand(selectedNeuronRand,1:numOdorsToPlot);
classResponsesPt2 = pirClassRand(selectedNeuronRand,1:numOdorsToPlot);


%% SAVE EVERYTHING

save('part2ForImportanceFig_pInv.mat',...
    'NoC6',...
    'NoM6',...
    'NoR6',...
    'importanceLoc',...
    'importanceLocRand',...
    'importanceBound1a',...
    'importanceBound1b',...
    'importanceBound1aRand',...
    'importanceBound1bRand',...
    'importanceSparseBound',...
    'classFourthOrderFull',...
    'classFourthOrderFullRand',...
    'classFourthOrderRand',...
    'classFourthOrderRandRand',...
    'classFourthOrderNoBottom',...
    'classFourthOrderNoBottomRand',...
    'classFourthOrderNoTop',...
    'classFourthOrderNoTopRand',...
    'selectedNeuronRand',...
    'randResponsesPt2',...
    'classResponsesPt2',...
    'numOdorsToPlot');








