
function calculateFig1parts
global Ny Sp Sx figureOne classLevel1 classLevel2 class1Color class2Color

% split calculation from plotting in makeFig1.m on 12/25/17

initializeGlobals_Server
tmpTime = clock;
rng(1e5*tmpTime(end)); %seed random number generator

No = 100;
NoTot = 3*No;
xsSteps = [.25, .5, .75, 1];
xsMat = repmat(reshape((ones(No/4,1)*xsSteps),1,No),Nx,3);
xsMat = xsMat>rand(size(xsMat));

Sclass = classLevel2*(ones(2*No)-eye(2*No))+eye(2*No);
Sclass(No+1:end,No+1:end) = classLevel1*(ones(No)-eye(No))+eye(No);
S = eye(NoTot);
S(No+1:end,No+1:end) = Sclass;

Sx = 0.2;
[glomRepTemp,odorSimilarity] = makeOdors('correlated',S,NoTot);
glomerulusRepRand = glomRepTemp(:,1:No).*xsMat(:,1:No);
similarityRand = -eps*(odorSimilarity(1:No)<0); 
glomerulusRepCorr2 = glomRepTemp(:,No+1:2*No).*xsMat(:,No+1:2*No);
similarityClass2 = -2 + odorSimilarity(No+1:2*No);
glomerulusRepCorr1 = glomRepTemp(:,2*No+1:end).*xsMat(:,2*No+1:end);
similarityClass1 = -2 + odorSimilarity(2*No+1:end);

[piriformRepRand,~,~,J] = makePiriform(glomerulusRepRand);
piriformRepClass1 = makePiriform(glomerulusRepCorr1,Sp,J);
piriformRepClass2 = makePiriform(glomerulusRepCorr2,Sp,J);

temp1 = 1/Ny*( double(piriformRepRand~=0)'*double(piriformRepRand~=0) );
pirRepFreqRandObserved = temp1(find(triu(ones(No))-eye(No)));
temp2 = 1/Ny^2*( sum(piriformRepRand~=0)'*sum(piriformRepRand~=0) );
pirRepFreqRandExpected = temp2(find(triu(ones(No))-eye(No)));

temp1 = 1/Ny*( double(piriformRepClass1~=0)'*double(piriformRepClass1~=0) );
pirRepFreqClass1Observed = temp1(find(triu(ones(No))-eye(No)));
temp2 = 1/Ny^2*( sum(piriformRepClass1~=0)'*sum(piriformRepClass1~=0) );
pirRepFreqClass1Expected = temp2(find(triu(ones(No))-eye(No)));

temp1 = 1/Ny*( double(piriformRepClass2~=0)'*double(piriformRepClass2~=0) );
pirRepFreqClass2Observed = temp1(find(triu(ones(No))-eye(No)));
temp2 = 1/Ny^2*( sum(piriformRepClass2~=0)'*sum(piriformRepClass2~=0) );
pirRepFreqClass2Expected = temp2(find(triu(ones(No))-eye(No)));

save('partsForFig1odorPairComparison.mat')

