
function plotModelCorrVsInput

initializeGlobals_SmallNy;

% this is the file that makes the model part of figure 2
% ------------------------
runRandPart = false;
runClassPart = false;
runMixPart = false;
runCorrPart = true;
% ------------------------


numSamples = 500; %8; % number of "datasets" per param value


numPts = 25; %50; %100; %number of points in each dataset
%p2Mat = NaN*ones( numSamples, length(f) ); % ratio of successful/total draws

%%
if runRandPart

    corrMatRand = NaN*ones( numSamples, 1 );
    
    parfor k=1:numSamples
        glomRep = makeOdors('classWithoutTemplate',0,2);
        piriformRep = makePiriform(glomRep);
        isBothOdors = find( sum( piriformRep>0, 2)==2 );
        n = 0;
        cutpt = min(length(isBothOdors),numPts);
        sample = zeros(numPts,2);
        sample(1:cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
        if length(isBothOdors)<numPts
            display( 'Not Enough Points' );
            n = cutpt;
            while n<numPts
                glomRep = makeOdors('classWithoutTemplate',0,2);
                piriformRep = makePiriform(glomRep);
                isBothOdors = find( sum( piriformRep>0, 2)==2 );
                cutpt = min(length(isBothOdors),numPts-n);
                sample(n+1:n+cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
                n=n+cutpt;
            end
            
        end
        
        C = corrcoef( sample(:,1), sample(:,2) );
        corrMatRand(k) = C(2);
    end
    

    db = 0.05;
    cBins = -1:db:1;
    h = hist(corrMatRand,cBins)/length(corrMatRand)/db;
    h1 = figure; bar(cBins,h, 'k' )
    box off
    set(gcf,'color','w')
    pos=get(gcf,'Position');
    pos(4)=pos(4)/2;
    set(gcf,'Position',pos)
    xlim([-1 1])
    ylabel('Probability Density','fontsize',12);
    xlabel('Correlation Coefficient','fontsize',12);
    saveas(h1,'~/Dropbox/_axellab/_Dan/figures/model/modelFig2 parts/modelRandDist.eps')
    save('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatmix.mat','corrMatRand','cBins','db');
end

%%
if runClassPart
    f = 0:.05:1;
    corrMat = NaN*ones( numSamples, length(f) );

    for j=1:length(f)
        display([num2str(j),' of ',num2str(length(f))])
        if f(j)==0.1
            if exist('corrMatRand','var')
                corrMat(:,j) = corrMatRand;
                continue;
            end
        end
        parfor k=1:numSamples
            glomRep = makeOdors('classWithoutTemplate',f(j),2);
            piriformRep = makePiriform(glomRep);
            isBothOdors = find( sum( piriformRep>0, 2)==2 );
            n = 0;
            cutpt = min(length(isBothOdors),numPts);
            sample = zeros(numPts,2);
            sample(1:cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
            if length(isBothOdors)<numPts
                display( 'Not Enough Points' );
                n = cutpt;
                while n<numPts
                    glomRep = makeOdors('classWithoutTemplate',f(j),2);
                    piriformRep = makePiriform(glomRep);
                    isBothOdors = find( sum( piriformRep>0, 2)==2 );
                    cutpt = min(length(isBothOdors),numPts-n);
                    sample(n+1:n+cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
                    n=n+cutpt;
                end
                
            end
            
            C = corrcoef( sample(:,1), sample(:,2) );
            corrMat(k,j) = C(2);
        end
        
    end
    
    h1 = figure; errorbar(f, mean(corrMat) , std(corrMat) , 'g.' )
    box off
    xlabel('Fraction Input Similarity','fontsize',12);
    ylabel('Correlation Coefficient','fontsize',12);
    saveas(h1,'~/Dropbox/_axellab/_Dan/figures/model/modelFig2 parts/modelDist.eps')
    save('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatClass.mat','corrMat','f');
    
    %h2 = figure; errorbar(f, mean(p2Mat) , std(p2Mat) , '.' )
    %xlabel('fraction input overlap','fontsize',12);
    %ylabel('P_{2}','fontsize',12);
    
end


%%
if runMixPart
    f = 0:.05:1;
    corrMat = NaN*ones( numSamples, length(f) );

    for j=1:length(f)
        display([num2str(j),' of ',num2str(length(f))])
        
        parfor k=1:numSamples
            glomRep = makeOdors('singleMix',f(j),2);
            piriformRep = makePiriform(glomRep);
            isBothOdors = find( sum( piriformRep>0, 2)==2 );
            n = 0;
            cutpt = min(length(isBothOdors),numPts);
            sample = zeros(numPts,2);
            sample(1:cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
            if length(isBothOdors)<numPts
                display( 'Not Enough Points' );
                n = cutpt;
                while n<numPts
                    glomRep = makeOdors('singleMix',f(j),2);
                    piriformRep = makePiriform(glomRep);
                    isBothOdors = find( sum( piriformRep>0, 2)==2 );
                    cutpt = min(length(isBothOdors),numPts-n);
                    sample(n+1:n+cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
                    n=n+cutpt;
                end
                
            end
            
            C = corrcoef( sample(:,1), sample(:,2) );
            corrMat(k,j) = C(2);
        end
        
    end
    
    figure(h1); hold all
    errorbar(f, mean(corrMat) , std(corrMat) , 'b.' )
    xlabel('fraction input overlap','fontsize',12);
    ylabel('Correlation Coefficient','fontsize',12);
    saveas(h1,'~/Dropbox/_axellab/_Dan/figures/model/modelFig2 parts/modelDist.eps')
    save('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatMix.mat','corrMat','f');
    
    %h2 = figure; errorbar(f, mean(p2Mat) , std(p2Mat) , '.' )
    %xlabel('fraction input overlap','fontsize',12);
    %ylabel('P_{2}','fontsize',12);
    
end

%%
if runCorrPart
    f = 0:.05:1;
    corrMat = NaN*ones( numSamples, length(f) );
    fracSignificant = zeros(1, length(f));

    for j=1:length(f)
        display([num2str(j),' of ',num2str(length(f))])
        
        for k=1:numSamples
        %parfor k=1:numSamples
            glomRep = makeOdors('correlated',f(j),2);
            piriformRep = makePiriform(glomRep);
            isBothOdors = find( sum( piriformRep>0, 2)==2 );
            n = 0;
            cutpt = min(length(isBothOdors),numPts);
            sample = zeros(numPts,2);
            sample(1:cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
            if length(isBothOdors)<numPts
                display( 'Not Enough Points' );
                n = cutpt;
                while n<numPts
                    glomRep = makeOdors('correlated',f(j),2);
                    piriformRep = makePiriform(glomRep);
                    isBothOdors = find( sum( piriformRep>0, 2)==2 );
                    cutpt = min(length(isBothOdors),numPts-n);
                    sample(n+1:n+cutpt,:) = piriformRep( isBothOdors(1:cutpt), : );
                    n=n+cutpt;
                end
            else
                fracSignificant(j) = fracSignificant(j)+1;
            end
            
            C = corrcoef( sample(:,1), sample(:,2) );
            corrMat(k,j) = C(2);
        end
        
    end
    
    fracSignificant = fracSignificant/numSamples;
    
    %figure(h1); hold all
    %errorbar(f, mean(corrMat) , std(corrMat) , 'b.' )
    %xlabel('fraction input overlap','fontsize',12);
    %ylabel('Correlation Coefficient','fontsize',12);
    %saveas(h1,'~/Dropbox/_axellab/_Dan/figures/model/modelFig2 parts/modelDist.eps')
    save('~/Dropbox/_axellab/matlab/_dan/_modelDist/_modelCorr/corrMatCorr.mat','corrMat','f','fracSignificant');
    
    %h2 = figure; errorbar(f, mean(p2Mat) , std(p2Mat) , '.' )
    %xlabel('fraction input overlap','fontsize',12);
    %ylabel('P_{2}','fontsize',12);
    
end
