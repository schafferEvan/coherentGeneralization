
function [y,th,h,J,maxy] = makePiriform(x,pirSparseness,J,maxy,computeJ)
global Nx Ny Sc ScInh exInhBalance glomExpMean glomExpMax Sp NoR NoC NoM classOdorOverlap N4 S4 S4c templateOdor Jmax glomActMu glomActSig thMean thVar


% rearranged to reduce RAM pressure on 1/2/2017


%initializeGlobals;
th = thMean; %7.8;%9.8;%9.4; %7.3;

if ~exist('pirSparseness','var')
    pirSparseness=Sp;
elseif isempty(pirSparseness)
    pirSparseness=Sp;
end

makeJ = false;
if ~exist('J','var')
    makeJ = true;
elseif isempty(J)
    makeJ = true;
end

if makeJ
    %Cp = Sc*Nx; %glomerulus connections per piriform neuron
    
    %random version of Jp
    %JpTemp = (rand(Ny,Nx)<Sc);%/Ny^2*; %randn(Ny,Nx).*
    
    %fixed version of Jp
    [~,loc] = sort(rand(Ny,Nx),2,'descend');
    JpTemp = zeros(Ny,Nx);
    for j=1:Ny
        JpTemp(j,loc(j,1:Sc*Nx),:) = 1;
    end
    
    
    JnTemp = logical(zeros(Ny,Nx));
    %old version of Jn
    %pLocTemp = find(~JpTemp); % find candidate inhib synapses
    %[~,orderTemp] = sort(rand(size(pLocTemp)),'descend');
    %JnTemp(pLocTemp(orderTemp(1:Nx*Ny*ScInh))) = 1; %location of inhib synapses
    %new version of Jn
    
    %pLocTemp = find(~JpTemp); % find candidate inhib synapses
    pLocTemp = ~JpTemp(:);
    
    % random version of Jn
    pNew = Ny*Nx*ScInh/sum(pLocTemp); %adjusted probability of inh
    isCon = rand(length(pLocTemp),1) < pNew;
    JnTemp( pLocTemp & isCon ) = 1;%/Ny^2; %location of inhib synapses
    % JnTemp(pLocTemp(isCon)) = 1;%/Ny^2; %location of inhib synapses
    
    

    
    % old version (requires more ram, 1/4/17)
    %pNew = Ny*Nx*ScInh/length(pLocTemp); %adjusted probability of inh
    %isCon = rand(size(pLocTemp)) < pNew;
    %JnTemp(pLocTemp(isCon)) = 1;%/Ny^2; %location of inhib synapses
    
    % fixed version of Jn
    %   for j=1:Ny
    %   pLocTemp = find(~JpTemp(j,:)); % find candidate inhib synapses
    %   [~,orderTemp] = sort(rand(size(pLocTemp)),'descend');
    %   JnTemp(j,pLocTemp(orderTemp(1:Nx*ScInh))) = 1; %location of inhib synapses
    %   end
    
    %     pValTemp = lognrnd(0,1,Ny,Nx);
    %     while max(max(pValTemp))>Jmax
    %         pValTemp(pValTemp>Jmax) = lognrnd(0,1,size(pValTemp(pValTemp>Jmax)));
    %     end
    %
    %     nValTemp = lognrnd(0,1,Ny,Nx);
    %     while max(max(nValTemp))>Jmax
    %         nValTemp(nValTemp>Jmax) = lognrnd(0,1,size(nValTemp(nValTemp>Jmax)));
    %     end
    %     nValTemp = nValTemp*exInhBalance*(Sc/ScInh); %adjust synaptic strengths to balance
    %
    %     J = JpTemp.*pValTemp - JnTemp.*nValTemp;
    if nargin>=5
        if computeJ
            J = JpTemp - JnTemp*exInhBalance*(Sc/ScInh);
            clear JpTemp JnTemp
        end
    else
        J = JpTemp - JnTemp*exInhBalance*(Sc/ScInh);
        clear JpTemp JnTemp
    end
end
clear pLocTemp

if exist('J','var') && ~isempty(J)
    y = J*x; %sparse(CPconn*glomerulusRep); %unthresholded activation of piriform
else
    J = [];
    y = JpTemp*x - JnTemp*x*exInhBalance*(Sc/ScInh);
    clear JpTemp JnTemp
end
h = NaN;
sz = size(y);
No = sz(2);

if Ny<=10000
    temp = sort(reshape(y,sz(1)*sz(2),1),'descend');
    temp2 = zeros(size(temp));
    temp2(1)=temp(1);
    for j=2:length(temp2)
        temp2(j)=temp2(j-1)+temp(j);
    end
    temp2(temp<0)=inf; %prevents following line from choosing -th instead of th
    [~,loc] = min(abs(temp2-temp2(round(Ny*pirSparseness*No))/.95)); % --> 10% cutoff accounts for 95% of total
    thCheck = temp(loc); %sum of first Sx fraction of responses / .95 = target
    display(['Threshold = ',num2str(th), ' = ', num2str( th/exp(glomActMu+glomActSig^2/2) ), ' avg glom inputs. Exact value for ',num2str(100*Sp),'% = ',num2str(thCheck) ] )
end

if ~exist('maxy','var') || isempty(maxy)
    %maxy = max(max(y));
    
    % instead of normalizing to largest value, normalize to 99th percentile
    % for greater stability of results
    tempsort = sort(y,'descend');
    %xyMax = full(max(max(mouse2data),max(mouse1data)));
    maxy = full(tempsort( round(.01*length(tempsort)) ));
end

%thMat = thMean*ones(size(y)) + sqrt(thVar)*repmat( randn(Ny,1), 1, No );
%y = y-thMat;
y = y-th; % y = sign(y-th);
y(y<0)=0;
%y(y<th) = 0;
th = th/maxy; % SO "TH" STILL MEANS SOMETHING ******
y = sparse(y/maxy);  








% ----- OLD ---------------------------------------------------------------



% % % amp dependence version --------------------------------
% sz = size(x);
% No = sz(2);
% [~,connLoc] = sort(rand(Ny,Nx),2,'descend');
% J = connLoc <= Sc*(Nx);
% %J = (rand(Ny,Nx)<Sc);
% 
% SpLevel = .95; %.95;
% h = J*x; %sparse(CPconn*glomerulusRep); %unthresholded activation of piriform
% y = h + .2*randn(size(h));
% [~,loc] = sort(h,'descend');
% for j=1:No
%     y(loc(Ny*SpLevel+1:end,j),j) = 0; %10;
% end
% 


% % older --------------------------------------------
% J = rand(Ny,Nx)<Sc; %randn(Ny,Nx).*
% J = J-mean(mean(J));
% %J = randn(Ny,Nx).*(rand(Ny,Nx)<Sc);
% %J = J - mean(mean(J));
% 
% h = J*x; %sparse(CPconn*glomerulusRep); %unthresholded activation of piriform
% y = h.*(h>0);
% y = y/max(max(y));

% ORIGINAL ------------------------------------------
% temp = rand(Nx,Ny);
% [~,temp2] = sort(temp,'descend');
% idx2 = reshape(temp2(1:Cp,:),Cp*Ny,1); 
% idx1 = reshape(repmat(1:Ny,Cp,1),Cp*Ny,1);
% CPconn = sparse(idx1,idx2,1,Ny,Nx); %connection matrix from glomeruli to piriform
% piriformInput = sparse(CPconn*x); %unthresholded activation of piriform
% [temp1,temp2] = sort(piriformInput,'descend');
% idx2 = reshape(temp2(1:Ny*pirSparseness,:),Ny*(NoR+NoC+NoM)*pirSparseness,1); 
% valp = reshape(temp1(1:Ny*pirSparseness,:),Ny*(NoR+NoC+NoM)*pirSparseness,1); 
% idx1 = reshape(repmat(1:(NoR+NoC+NoM),Ny*pirSparseness,1),Ny*(NoR+NoC+NoM)*pirSparseness,1);
% y = sparse(idx2,idx1,valp,Ny,(NoR+NoC+NoM)); %thresholded representation in piriform




% OLD ---------------------------------------------------------

% % piriform (mouse 1)
% 
% temp = rand(Nx,Ny);
% [~,temp2] = sort(temp,'descend');
% idx2 = reshape(temp2(1:Cp,:),Cp*Ny,1); 
% idx1 = reshape(repmat(1:Ny,Cp,1),Cp*Ny,1);
% CPconn = sparse(idx1,idx2,1,Ny,Nx); %connection matrix from glomeruli to piriform
% clear temp temp2
% 
% % piriform (mouse 2)
% 
% temp = rand(Nx,Ny);
% [~,temp2] = sort(temp,'descend');
% idx2 = reshape(temp2(1:Cp,:),Cp*Ny,1); 
% idx1 = reshape(repmat(1:Ny,Cp,1),Cp*Ny,1);
% CPconn2 = sparse(idx1,idx2,1,Ny,Nx); %connection matrix from glomeruli to piriform
% clear temp temp2

% % glomeruli -> piriform (mouse 1)
% piriformInput = sparse(CPconn*glomerulusRep); %unthresholded activation of piriform
% [temp1,temp2] = sort(piriformInput,'descend');
% idx2 = reshape(temp2(1:Ny*pirSparseness,:),Ny*(NoR+NoC+NoM)*pirSparseness,1); 
% valp = reshape(temp1(1:Ny*pirSparseness,:),Ny*(NoR+NoC+NoM)*pirSparseness,1); 
% idx1 = reshape(repmat(1:(NoR+NoC+NoM),Ny*pirSparseness,1),Ny*(NoR+NoC+NoM)*pirSparseness,1);
% %display('THIS IS THE LINE I JUST CHANGED')
% piriformRep = sparse(idx2,idx1,valp,Ny,(NoR+NoC+NoM)); %thresholded representation in piriform
% %piriformRep = piriformInput;
% clear temp1 temp2 valp