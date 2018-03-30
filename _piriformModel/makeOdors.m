
function [x,similarity,overlap] = makeOdors(kind,value,No)
global Nx Sx glomExpMean glomExpMin glomExpMax templateOdor glomActMu glomActSig
%initializeGlobals;

similarity = ones(1,No);
overlap = nan(1,No);
if strcmp(kind,'classWithoutTemplate')
    fracClassGlom = value;
    cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
    %pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( round((1-fracClassGlom*Sx)*Nx), No)];
    pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( Nx-round(fracClassGlom*Sx*Nx), No)];
    x = rand(Nx,No) < pMat;
    x = makeXmagnitudes(x,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
    similarity = value*ones(1,No);
    
elseif strcmp(kind,'classWithTemplate')
    fracClassGlom = value;
    cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
    pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( round((1-fracClassGlom*Sx)*Nx), No)];
    x = rand(Nx,No) < pMat;
    x = makeXmagnitudes(x,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
    similarity = -ones(1,No);
    x(:,1) = templateOdor;
    
elseif strcmp(kind,'mix')
    % this gives a spectrum of mixes with *DIFFERENT* weights
    x = zeros(Nx,No);
    fracClassGlom = 0;
    cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
    pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( round((1-fracClassGlom*Sx)*Nx), No)];
    xtmp = rand(Nx,No) < pMat;
    xtmp = makeXmagnitudes(xtmp,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
    for j=0:No-1
        x(:,j+1) = sqrt(1-j/No)*xtmp(:,j+1) + sqrt(j/No)*templateOdor;
        similarity(j+1) = No-j; %j+1;
    end
    
elseif strcmp(kind,'singleMix')
    % this gives a set of mixes with the *SAME* weights
    x = zeros(Nx,No);
    fracClassGlom = 0;
    cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
    pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( round((1-fracClassGlom*Sx)*Nx), No)];
    xtmp = rand(Nx,No) < pMat;
    xtmp = makeXmagnitudes(xtmp,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
    for j=0:No-1
        x(:,j+1) = sqrt(1-value)*xtmp(:,j+1) + sqrt(value)*templateOdor;
        similarity(j+1) = value; %j+1;
    end
    
elseif strcmp(kind,'correlated')
    %pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( round((1-fracClassGlom*Sx)*Nx), No)];
    if length(value)==1
        fracClassGlom = value;
        cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
        pMat = [ones( round(fracClassGlom*Sx*Nx), No); cVal*ones( Nx-round(fracClassGlom*Sx*Nx), No)];
        similarity = value*ones(1,No);
    else
        pMat = nan(Nx,No);
        similarity = nan(1,No);
        for j=1:length(value)
            if j==1
                fracClassGlom = value(2,1);
            else
                fracClassGlom = value(j-1,j);
            end
            cVal = (1-fracClassGlom)*Sx/(1 - fracClassGlom*Sx);
            pMat(:,j) = [ones( round(fracClassGlom*Sx*Nx), 1); cVal*ones( Nx-round(fracClassGlom*Sx*Nx), 1)];
            similarity(j) = fracClassGlom;
        end
    end
    xtmp = rand(Nx,No) < pMat;
    x = makeXmagnitudesCorr(xtmp,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig,value);
end



    



%% OLD

% 
% elseif strcmp(kind,'variableOverlapOLD')
%     fracClassGlom = value;
%     numClassGlom = round(fracClassGlom*Nx*Sx);
%     numNonClassGlom = round(Nx*Sx - numClassGlom);
% 
%     x = zeros(Nx,No);
%     temp = rand(Nx,1);
%     [~,loc] = sort(temp,'descend');
%     x(loc(1:Sx*Nx),1) = 1;         %this one is the template
%     templateLocs = find(x(:,1));
%     notTemplateLocs = find(~x(:,1));
%     
%     for j=2:No
%         tempClass = rand(Nx*Sx,1);
%         tempNoClass = rand(Nx*(1-Sx),1);
%         [~,locClass] = sort(tempClass,'descend');
%         [~,locNoClass] = sort(tempNoClass,'descend');
%         x(templateLocs(locClass(1:numClassGlom)),j) = 1;
%         x(notTemplateLocs(locNoClass(1:numNonClassGlom)),j) = 1;
%     end
%     x = makeXmagnitudes(x,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
%     similarity = value*ones(1,No);
% 
% elseif strcmp(kind,'fig1aOLD')
%     
%     %x0 = rand(Nx,1);
%     x = zeros(Nx,No);
%     x1 = rand(Nx,No);
%     for j=0:No-1
%         
%         tmp = x1(:,j+1);
%         similarity(j+1) = -(j+1);
%         
%         [~,loc] = sort(tmp,'descend');
%         x(loc(1:Sx*Nx),j+1) = 1;
%     end
%     x = makeXmagnitudes(x,glomExpMin,glomExpMean,glomExpMax,No,Nx,glomActMu,glomActSig);
%     
    
% if strcmp(kind,'ampDependence')
%     fracPublicGlom = value;
% 
%     x = zeros(Nx,No);
%     temp = rand(Nx,1);
%     [~,locA] = sort(temp,'descend');
%     x(locA(1:Sx*Nx),1) = 1;
%     for j=2:No
%         temp = rand(Nx,1);
%         if fracPublicGlom>=0 % if positive (less than chance overlap), remove some options from x1's glomeruli
%             temp(locA(1:round(fracPublicGlom*Sx*Nx))) = 0;
%         else                 % if negative (greater than chance overlap), remove some options from NOT x1's glomeruli
%             %temp(locA(end-round(1-fracPublicGlom*(1-Sx)*Nx):end)) = 0;
%             temp(locA(end-round(-1-fracPublicGlom*(1-Sx)*Nx):end)) = 0;
%         end
%         [~,loc] = sort(temp,'descend');
%         x(loc(1:Sx*Nx),j) = 1;
%     end
%     
%     
%     x = x*glomExpMin + x.*random('exp',glomExpMean,Nx,No);
%     while sum(sum(x>glomExpMax))>0
%         x(x>glomExpMax)=random('exp',glomExpMean,size(x(x>glomExpMax))); %0;
%     end
%     x = x/glomExpMax;




% ODORS ---------------------------------------------------------

% %odors outside class
% temp = rand(NoR,Nx);%constraint is intentionally flipped
% [~,temp2] = sort(temp,2,'descend'); %constraint is intentionally flipped
% idx2 = reshape(temp2(:,1:Nx*Sx)',round(Nx*NoR*Sx),1); %constraint is intentionally flipped
% idx1 = reshape(repmat(1:NoR,Nx*Sx,1),round(Nx*NoR*Sx),1);%constraint is intentionally flipped
% glomerulusMap = sparse(idx1,idx2,1,NoR,Nx); %glomeruli activated for each odor
% glomerulusMap = glomerulusMap';
% clear temp temp2
% 
% % odors in class
% 
% NxOverlap = Nx*Sx*classOdorOverlap;
% NxNonOverlap = Nx - Nx*Sx*classOdorOverlap;
% temp = rand(NoC,NxNonOverlap);%constraint is intentionally flipped
% [~,temp2] = sort(temp,2,'descend'); %constraint is intentionally flipped
% idx2 = reshape(temp2(:,1:(Nx*Sx-NxOverlap))',(Nx*Sx-NxOverlap)*NoC,1); %constraint is intentionally flipped
% idx1 = reshape(repmat(1:NoC,(Nx*Sx-NxOverlap),1),(Nx*Sx-NxOverlap)*NoC,1);%constraint is intentionally flipped
% classGlomerulusMap = sparse(idx1,idx2,1,NoC,NxNonOverlap); %glomeruli activated for each odor
% classGlomerulusMap = classGlomerulusMap';
% classGlomerulusMap = [ones(NxOverlap,NoC); classGlomerulusMap];
% clear temp temp2
% 
% glomerulusMap = [classGlomerulusMap, glomerulusMap];

% % odors -> glomeruli
% 
% temp = random('exp',glomExpMean,Nx,NoR+NoC); %NOTE: this is the only place that shouldn't have NoM (mixtures defined below)
% %temp = ones(size(random('exp',glomExpMean,Nx,NoR+NoC))); %NOTE: this is the only place that shouldn't have NoM (mixtures defined below)
% 
% glomRepComplete = false;
% while ~glomRepComplete
%     remainingEntries = sum(sum(temp>glomExpMax));
%     if ~remainingEntries
%         glomRepComplete = true;
%     else
%         %temp(temp>glomExpMax) = ones(size(random('exp',glomExpMean,remainingEntries,1)));
%         temp(temp>glomExpMax) = random('exp',glomExpMean,remainingEntries,1);
%     end
% end
% %temp(temp>glomExpMax) = 0;
% glomerulusRep = sparse(temp.*glomerulusMap); %representation in glomeruli
% 
% % add mixture odors
% mixtemp = repmat(1/NoM:1/NoM:1,1000,1).*repmat(glomerulusRep(:,trainedOdorNum),1,NoM)...
%     + repmat(1:-1/NoM:1/NoM,1000,1).*glomerulusRep(:,end-NoM+1:end);
% %mixtemp = abs(repmat(1/NoM:1/NoM:1,1000,1).*repmat(glomerulusRep(:,trainedOdorNum),1,NoM)...
% %    + repmat(1:-1/NoM:1/NoM,1000,1).*repmat(glomerulusRep(:,trainedOdorNum),1,NoM).*randn(size(repmat(glomerulusRep(:,trainedOdorNum),1,NoM))));
% glomerulusRep = [glomerulusRep,mixtemp];
% clear mixtemp