

global Nx Ny Sx Sc ScInh exInhBalance glomExpMean glomExpMin glomExpMax Sp NoR NoC NoM classOdorOverlap N4 S4 S4c templateOdor Jmax glomActMu glomActSig thMean thVar outlierThresh classLevel1 classLevel2 greenColorAdjust

% PARAMETERS ---------------------------------------------------------
randn(1000,1); %for random number generator

Nx = 1000; %number of glomeruli
Ny = 100000; %number of piriform neurons  NOTE: THIS IS CORRECT NUMBER OF INPUTS
N4 = 50;

Sx = 0.1; % statement at right isn't true b/c of glomExpMin %.14; %fraction of glomeruli activated by an odor -- 14% exp gives 10% accounting for 95% of activity
Sc = .2;%.1; %.2; %.1; %.2; %prob of conn from glom to pir
ScInh = .4;%0.2;
exInhBalance = 1; % <1 favors Ex, >1 favors Inh.  NOTE: THIS IS DEPRECATED.  KEEP THIS FIXED AND IGNORE IT


glomActMu = .1;%0.07;%0.1; %"mu" of lognormal dist
glomActSig = .5;%0.7;%0.5;%"sig" of lognormal dist
glomExpMean = 1;    % DEPRECATED (TO BE REPLACED BY GLOMACTMU)
glomExpMin = .5;    % DEPRECATED (TO BE REPLACED BY GLOMACTSIG)
glomExpMax = 20; %5*glomExpMean;

Sp = 0.062; %.1; %.08; %.1;
S4 = Sp;
S4c = .1;
Jmax = 20; %4;

thMean = 11.87;%11.5;%/Ny^2;%8.4; %piriform threshold
thVar = 0;%1.8;

NoR = 200; %number of test odors
NoC = 200; %number of odors in class
NoM = 2000; %number of odors in mix
classOdorOverlap = 1; %0.98; %fraction of shared Glomeruli in a class

classLevel1 = 0.7;
classLevel2 = 0.3;
greenColorAdjust = 0.25;


%outlierThresh = inf; %threshold above which to consider datapoints outliers for figure 2.
outlierThresh = inf; %1-1e-15; %1-1e-14; %

% -- template -----
if ~exist('templateOdor','var') || isempty(templateOdor)
    templateOdor = makeOdors('classWithoutTemplate',classLevel1,1);
    %     templateOdor = zeros(Nx,1);
    %     templateOdor(1:Nx*Sx*classOdorOverlap) = 1; %random('exp',glomExpMean,Nx*Sx*classOdorOverlap,1);
    %     templateOdor = makeXmagnitudes(templateOdor,glomExpMin,glomExpMean,glomExpMax,1,Nx,glomActMu,glomActSig);
    %     %templateOdor(templateOdor>glomExpMax)=random('exp',glomExpMean,size(templateOdor(templateOdor>glomExpMax)));
    %     %templateOdor = templateOdor/glomExpMax;
end
%------------------
display('initialization complete')
