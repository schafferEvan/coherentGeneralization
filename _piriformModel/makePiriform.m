
function [y,th,h,J,maxy,thNoNorm] = makePiriform(x,pirSparseness,J,maxy,computeJ,th)
global Nx Ny Sc ScInh exInhBalance Sp glomActMu glomActSig thMean


if nargin<6
    th = thMean;
    useExactTh = true;
else
    useExactTh = false;
end
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
    JpTemp = (rand(Ny,Nx)<Sc);
    JnTemp = logical(zeros(Ny,Nx));
    pLocTemp = ~JpTemp(:);
    
    pNew = Ny*Nx*ScInh/sum(pLocTemp); %adjusted probability of inh
    isCon = rand(length(pLocTemp),1) < pNew;
    JnTemp( pLocTemp & isCon ) = 1;%/Ny^2; %location of inhib synapses
    
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
    y = J*x;
else
    J = [];
    y = JpTemp*x - JnTemp*x*exInhBalance*(Sc/ScInh);
    clear JpTemp JnTemp
end
h = NaN;
sz = size(y);
No = sz(2);

if Ny<=10000 && useExactTh
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
    
    if Ny<=10000
        th = thCheck;
        display('using exact theta for small Ny')
    end
    
end

if ~exist('maxy','var') || isempty(maxy)
    
    % instead of normalizing to largest value, normalize to 99th percentile
    % for greater stability of results
    tempsort = sort(y,'descend');
    maxy = full(tempsort( round(.01*length(tempsort)) ));
end

y = y-th;
y(y<0)=0;
thNoNorm = th;
th = th/maxy;
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