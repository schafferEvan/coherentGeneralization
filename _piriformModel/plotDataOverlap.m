
function [nCellsA,nCellsB,nCellsTot,observedOverlap,expectedOverlap,compNames] = plotDataOverlap
% PLOTS DATA CORRELATION



global onlyAnalyzeCoactiveCells ignoreZeroBin wholeOrSampleDataset optimizeLog figureOne outlierThresh

%----------
onlyAnalyzeCoactiveCells = true;
ignoreZeroBin = true; % with adaptive binning, this has to be true
optimizeLog = false;
%----------




[responses,Ncells,Nodors,odorNames,isMix,fracResponsive] = danLoadBigDataset;
numPairs = nchoosek(Nodors,2);
bSparseness = NaN*ones(numPairs,1);
aSparseness = NaN*ones(numPairs,1);
nCellsBoth = NaN*ones(numPairs,1);
numOutliers = NaN*ones(numPairs,1);
observedOverlap = NaN*ones(numPairs,1);
expectedOverlap = NaN*ones(numPairs,1);
nCellsA = NaN*ones(numPairs,1);
nCellsB = NaN*ones(numPairs,1);
nCellsTot = NaN*ones(numPairs,1);
bothOdorsAboveThreshold = NaN*ones(numPairs,1);
pairIsSingleOdors = NaN*ones(numPairs,1);
dataOverlap = NaN*ones(numPairs,1);
compNames = cell(numPairs,1);
monoOdors = find(~isMix);
mixOdors = find(isMix);
componentIsShared = NaN*ones(numPairs,1);
singleOdorCompCount = 0;
sharedMixCompCount = 0;
unsharedMixCompCount = 0;
cellThresh = 10;

odorPartsList = struct('includes',{},'without',{});
for j=1:length(odorNames)
    strtstp = strfind(odorNames{j},'_');
    adds = strfind(odorNames{j},'+');
    subs = strfind(odorNames{j},'-'); % below assumes only special case of subs
    isdate = strfind(odorNames{j},'090115');
    if ~isempty(adds)
        if isempty(subs)
            odorPartsList(j).includes{1} = odorNames{j}(strtstp(1)+1:adds(1)-1);
            if isempty(isdate)
                for k=2:length(adds)
                    odorPartsList(j).includes{k} = odorNames{j}(adds(k-1)+1:adds(k)-1);
                end
                odorPartsList(j).includes{length(adds)+1} = odorNames{j}(adds(end)+1:strtstp(2)-1);
            end
        else
            odorPartsList(j).includes{1} = odorNames{j}(strtstp(1)+1:subs(1)-1);
            %odorPartsList(j).includes{2} = odorNames{j}(adds(1)+1:strtstp(2)-1);
            odorPartsList(j).without{1} = odorNames{j}(subs(1)+1:adds(1)-1);
        end
    else
        odorPartsList(j).includes = {odorNames{j}(strtstp(1)+1:strtstp(2)-1)};
    end
end

ctr = 0;
for i=1:length(monoOdors)
    oName1 = [];
    for ii=1:length(odorPartsList(monoOdors(i)).includes)
        oName1 = [oName1,odorPartsList(monoOdors(i)).includes{ii}];
        if ii<length(odorPartsList(monoOdors(i)).includes)
            oName1 = [oName1,'+'];
        end
        if ~isempty(odorPartsList(monoOdors(i)).without)
            oName1 = [oName1,' w/o ',odorPartsList(monoOdors(i)).without{1}];
        end
    end
    for j=i+1:length(monoOdors)
        ctr = ctr + 1;
        
        pairIsSingleOdors(ctr) = 1;
        oName2 = [];
        for ii=1:length(odorPartsList(monoOdors(j)).includes)
            oName2 = [oName2,odorPartsList(monoOdors(j)).includes{ii}];
            if ii<length(odorPartsList(monoOdors(j)).includes)
                oName2 = [oName2,'+'];
            end
            if ~isempty(odorPartsList(monoOdors(j)).without)
                oName2 = [oName2,' w/o ',odorPartsList(monoOdors(j)).without{1}];
            end
        end
        compNames{ctr} = [oName1,' vs ',oName2];
        [respA,respB,aSparseness(ctr),bSparseness(ctr)] = danGetBinnedOdorPair(monoOdors(i),monoOdors(j),responses,isMix);
        
        
        bothActive = (respB>0) & (respA>0);
        aActive = (respA>0);
        bActive = (respB>0);
        %isValid = sum( ~isnan(respA) & ~isnan(respB) );
        %if length(respA)>cellThresh
        if (sum(aActive)>cellThresh) && (sum(bActive)>cellThresh)
            singleOdorCompCount = singleOdorCompCount + 1;
        end
        
        %muA = mean(respA(bothActive));
        %muB = mean(respB(bothActive));
        %sigA = std(respA(bothActive));
        %sigB = std(respB(bothActive));
        aPts = respA(bothActive);
        bPts = respB(bothActive);
        mvMu = mean([aPts,bPts]);
        mvSig = cov(respA(bothActive),respB(bothActive));
        %isOutlier = (respA(bothActive)>muA+outlierTh*sigA) & ...
        %    (respB(bothActive)>muB+outlierTh*sigB);
        mahalD = zeros(size(aPts)); % mahalanobis distance
        isOutlier = zeros(size(aPts));
        for jj=1:length(mahalD)
            Y = [aPts(jj),bPts(jj)];
            mahalD(jj) = (Y-mvMu)/mvSig*(Y-mvMu)'; %(Y-mvMu)*inv(mvSig)*(Y-mvMu)';
            isOutlier(jj) = chi2cdf(mahalD(jj),2) > outlierThresh;
        end
        
        
        numOutliers(ctr) = sum(isOutlier);
        nCellsBoth(ctr) = sum(bothActive) - sum(isOutlier);
        nCellsA(ctr) = sum(aActive) - sum(isOutlier);
        nCellsB(ctr) = sum(bActive) - sum(isOutlier);
        nCellsTot(ctr) = length(respA) - sum(isOutlier);
        bothOdorsAboveThreshold(ctr) = (nCellsA(ctr)>cellThresh) && (nCellsB(ctr)>cellThresh);
        if nCellsBoth(ctr)>cellThresh
            if isempty(strfind(oName1,'uri')) && isempty(strfind(oName2,'uri'))
                observedOverlap(ctr) = nCellsBoth(ctr)/nCellsTot(ctr);
                expectedOverlap(ctr) = (nCellsA(ctr)/nCellsTot(ctr)) * (nCellsB(ctr)/nCellsTot(ctr));
            end
        end
        if strcmp( compNames{ctr}, 'pin vs eug')
            ex1.respA = respA;
            ex1.respB = respB;
            ex1.cc = dataOverlap(ctr);
        elseif strcmp( compNames{ctr}, 'hex vs oct')
            ex2.respA = respA;
            ex2.respB = respB;
            ex2.cc = dataOverlap(ctr);
        end
        
    end
end

firstMix = ctr+1;
display(firstMix)
for i=1:length(mixOdors)
    oName1 = [];
    for ii=1:length(odorPartsList(mixOdors(i)).includes)
        oName1 = [oName1,odorPartsList(mixOdors(i)).includes{ii}];
        if ii<length(odorPartsList(mixOdors(i)).includes)
            oName1 = [oName1,'+'];
        end
        if ~isempty(odorPartsList(mixOdors(i)).without)
            oName1 = [oName1,' w/o ',odorPartsList(mixOdors(i)).without{1}];
        end
    end
    
    for j=[1:mixOdors(1)-1,mixOdors(i)+1:Nodors] %note:skipping prev mixOdors prevents double-counting
        ctr = ctr + 1;

        
        oName2 = [];
        for ii=1:length(odorPartsList(j).includes)
            oName2 = [oName2,odorPartsList(j).includes{ii}];
            if ii<length(odorPartsList(j).includes)
                oName2 = [oName2,'+'];
            end
            if ~isempty(odorPartsList(j).without)
                oName2 = [oName2,' w/o ',odorPartsList(j).without{1}];
            end
        end
        
        compNames{ctr} = [oName1,' vs ',oName2];
        for k=1:length(odorPartsList(mixOdors(i)).includes)
            if componentIsShared(ctr) == 1
                break
            end
            for l=1:length(odorPartsList(j).includes)
                if strcmp( odorPartsList(mixOdors(i)).includes(k),...
                        odorPartsList(j).includes(l) )
                    componentIsShared(ctr) = 1;
                    break
                else
                    componentIsShared(ctr) = 0;
                end
            end
        end
        
        [respA,respB,aSparseness(ctr),bSparseness(ctr)] = danGetBinnedOdorPair(mixOdors(i),j,responses,isMix);
        
        bothActive = (respB>0) & (respA>0);
        aActive = (respA>0);
        bActive = (respB>0);
        
        if (sum(aActive)>cellThresh) && (sum(bActive)>cellThresh)
            if componentIsShared(ctr) == 1
                sharedMixCompCount = sharedMixCompCount + 1;
            else
                unsharedMixCompCount = unsharedMixCompCount + 1;
            end
        end
        
        %muA = mean(respA(bothActive));
        %muB = mean(respB(bothActive));
        %sigA = std(respA(bothActive));
        %sigB = std(respB(bothActive));
        aPts = respA(bothActive);
        bPts = respB(bothActive);
        mvMu = mean([aPts,bPts]);
        mvSig = cov(respA(bothActive),respB(bothActive));
        %isOutlier = (respA(bothActive)>muA+outlierTh*sigA) & ...
        %    (respB(bothActive)>muB+outlierTh*sigB);
        mahalD = zeros(size(aPts)); % mahalanobis distance
        isOutlier = zeros(size(aPts));
        for jj=1:length(mahalD)
            Y = [aPts(jj),bPts(jj)];
            mahalD(jj) = (Y-mvMu)/mvSig*(Y-mvMu)'; %(Y-mvMu)*inv(mvSig)*(Y-mvMu)';
            isOutlier(jj) = chi2cdf(mahalD(jj),2) > outlierThresh;
        end
        
        
        numOutliers(ctr) = sum(isOutlier);
        nCellsA(ctr) = sum(aActive) - sum(isOutlier);
        nCellsB(ctr) = sum(bActive) - sum(isOutlier);
        nCellsBoth(ctr) = sum(bothActive) - sum(isOutlier);
        if nCellsBoth(ctr)>cellThresh
            if isempty(strfind(oName1,'mix')) && isempty(strfind(oName2,'mix'))
                %C = corrcoef( respA(bothActive), respB(bothActive) );
                %dataOverlap(ctr) = C(2);
            end
        end
    end
end



figure(figureOne);
%subplot(3,12,[13:16,25:28]); hold on
axes('position',[.1  .08  .37  .37]); hold on


plot(expectedOverlap,observedOverlap,'.','markersize',10,'color',[0 0 0]);

box off;
xlabel('Chance Frequency','fontsize',12)
ylabel('Observed Frequency','fontsize',12)
%xlim([0 .02]);
%ylim([0 .02]);
set(gca,'XTick',0:.05:.2,'fontsize',12);
set(gca,'YTick',0:.02:.1,'fontsize',12);
lineTop = min( max(xlim),max(ylim) );
plot([0 lineTop],[0 lineTop],'color',[.5 .5 .5],'linewidth',3)

