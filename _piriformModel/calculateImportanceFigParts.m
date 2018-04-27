
%function makeImportanceFig
global Nx Ny Sx Sc glomExpMean glomExpMax Sp NoR NoC NoM classOdorOverlap N4 S4 S4c templateOdor paperFig

rand(1000);
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
% [glomRepTempMix,odorSimilarity] = makeOdors('mix',0,NoM6);
% mixtureOdorSimilarity = odorSimilarity;
% similarityMeasureFull = [similarityMeasureFull, -1+odorSimilarity/NoM6];
% glomerulusRep = [glomRepTempClass, glomRepTempRand];%, glomRepTempMix];

glomRepMouse2 = makeOdors('correlated',0,NoR6);

% PIRIFORM -------------------------------------------------------

[pirClass,~,~,J] = makePiriform(glomRepTempClass,Sp,[],maxy);  % mouse 1
pirRand = makePiriform(glomRepTempRand,Sp,J,maxy);  % mouse 1
pirClassRand = makePiriform(glomRepMouse2,Sp,J,maxy); % mouse 2


%%  CLASS DISCRIMINANT & IMPORTANCE (class example) -------------------------

%pirClass = piriformRep(:,1:NoC6);
%pirRand = piriformRep(:,NoC6+1:NoC6+NoR6);
mc = mean(pirClass,2);
mr = mean(pirRand,2);
%sc = full(sum(pirClass>0,2));
%sr = full(sum(pirRand>0,2));
%mcMat = repmat( mc , 1,200);
%mrMat = repmat( mr , 1,200);
%withinCov = (pirClass-mcMat)*(pirClass'-mcMat') + (pirRand-mrMat)*(pirRand'-mrMat');
%FisherLDweights = withinCov \ (mc-mr);
%FisherLDweights = FisherLDweights/abs(sum(FisherLDweights));

importance = (mc-mr);% ./ mr;% ./ (mc+mr);

[~,importanceLoc] = sort(importance,'descend');
importanceBound1 = round(.04*Ny);
importanceBound2 = find(importance(importanceLoc)>0, 1, 'last'); %Ny;
importanceSparseBound = S4c; %.99;
temp = rand(size(mc))<importanceSparseBound;

classWeightsFull = (mc-mr).*( (mc-mr)>0 );         % *** readout can't go negative.  old version was just (mc-mr) ***
classWeightsNoTop = classWeightsFull;
classWeightsNoTop(importanceLoc(1:importanceBound1)) = 0;
classWeightsNoBottom = classWeightsFull;
classWeightsNoBottom(importanceLoc(importanceBound1+1:importanceBound2)) = 0;
classWeightsRand = classWeightsFull;
classWeightsRand(temp) = 0;

classFourthOrderFull = full(classWeightsFull'*[pirClass,pirRand]);
classFourthOrderNoTop = full(classWeightsNoTop'*[pirClass,pirRand]);
classFourthOrderNoBottom = full(classWeightsNoBottom'*[pirClass,pirRand]);
classFourthOrderRand = full(classWeightsRand'*[pirClass,pirRand]);

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

selectedNeuron = importanceLoc(importanceBound1);
classResponsesPt1 = pirClass(selectedNeuron,1:numOdorsToPlot);
randResponsesPt1 = pirRand(selectedNeuron,1:numOdorsToPlot);

save('part1ForImportanceFig.mat','selectedNeuron','classResponsesPt1','randResponsesPt1');
clear pirClass
    
%% CLASS DISCRIMINANT & IMPORTANCE (RANDOM example) -------------------------

%pirClassRand = pirClassRand;%(:,NoC6+1:NoC6+NoR6);
%pirRand = piriformRep(:,NoC6+1:NoC6+NoR6);
mcr = mean(pirClassRand,2);
%mr = mean(pirRand,2);

importanceRand = (mcr-mr);% ./ mr;% ./ (mcr+mr);

[~,importanceLocRand] = sort(importanceRand,'descend');
importanceBound1Rand = round(.04*Ny);
importanceBound2Rand = find(importanceRand(importanceLocRand)>0, 1, 'last'); %Ny;
%importanceSparseBound = Sc; %.99;
temp = rand(size(mcr))<importanceSparseBound;

classWeightsFullRand = (mcr-mr).*( (mcr-mr)>0 );           % *** readout can't go negative.  old version was just (mcr-mr) ***
classWeightsNoTopRand = classWeightsFullRand;
classWeightsNoTopRand(importanceLocRand(1:importanceBound1Rand)) = 0;
classWeightsNoBottomRand = classWeightsFullRand;
classWeightsNoBottomRand(importanceLocRand(importanceBound1Rand+1:importanceBound2Rand)) = 0;
classWeightsRandRand = classWeightsFullRand;
classWeightsRandRand(temp) = 0;

%pirRepTemp = [pirClassRand,pirRand];
%pirRepTemp(:,1:NoC6) = pirClassRand;%(:,NoC6+1:NoC6+NoR6);
classFourthOrderFullRand = full(classWeightsFullRand'*[pirClassRand,pirRand]);
classFourthOrderNoTopRand = full(classWeightsNoTopRand'*[pirClassRand,pirRand]);
classFourthOrderNoBottomRand = full(classWeightsNoBottomRand'*[pirClassRand,pirRand]);
classFourthOrderRandRand = full(classWeightsRandRand'*[pirClassRand,pirRand]);

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

selectedNeuronRand = importanceLocRand(importanceBound1Rand);
randResponsesPt2 = pirRand(selectedNeuronRand,1:numOdorsToPlot);
classResponsesPt2 = pirClassRand(selectedNeuronRand,1:numOdorsToPlot);


%% SAVE EVERYTHING

save('part2ForImportanceFig.mat',...
    'NoC6',...
    'NoM6',...
    'NoR6',...
    'importanceLoc',...
    'importanceLocRand',...
    'importanceBound1',...
    'importanceBound1Rand',...
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








