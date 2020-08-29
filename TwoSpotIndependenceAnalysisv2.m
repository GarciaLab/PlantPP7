function TwoSpotIndependenceAnalysisv2(Prefix,DynamicsResultsPath,pickedFrame,fractionCompetent)


%% Load stuff
PrefixPath = [DynamicsResultsPath '/' Prefix];
SingleTracesPath = [DynamicsResultsPath '/' Prefix '/' ['ResultsFigures_' Prefix]];
load([PrefixPath '/' 'CompiledParticles.mat'],'CompiledParticles')
%in some versions of the code CompiledParticles can be a struct within a cell
if iscell(CompiledParticles) 
    CompiledParticles = CompiledParticles{1};
end
load([PrefixPath '/' 'Spots.mat'])
load([PrefixPath '/' 'Particles.mat'],'Particles')
load([PrefixPath '/' 'FrameInfo.mat'])
AbsTime = [FrameInfo.Time]./60; %in minutes
load([PrefixPath '/' 'FrameInfo.mat'])
load([PrefixPath '/' 'Ellipses.mat'])

% asign spots to nuclei from scratch based on consistent proximity
[NucleiPerFrame, CompiledParticles] = ...
    assignParticlesToEllipses(Ellipses,CompiledParticles,FrameInfo);
% find homolog allele spots
CompiledParticles = findSisterParticle(CompiledParticles);

%% Obtain p

% In each frame we count the number of spots and the number of nuclei to
% calculate the probability 'p' that an allele is actively transcribing given
% by: a/2N where a = number of spots; N = number of nuclei.

%get the number of nuclei per frame
FigFileName = 'CellsPerFrame.fig';
F1 = open([SingleTracesPath '/' FigFileName]);
YdataObjs = findobj(F1,'-property','YData');
N = YdataObjs(2).YData; %nuclei per frame

%get the number of spots per frame
FigFileName = 'SpotNumber.fig';
F2 = open([SingleTracesPath '/' FigFileName]);
YdataObjs = findobj(F2,'-property','YData');
a = YdataObjs(1).YData; %spots per frame
    
% calculate 'p', the probability that an allele is on
p = a./(2*N);

%close all
%% Obtain the expected number of nuclei with 0,1 and 2 spots in each frame

pZeroSpots = (1-p).^2;
pOneSpot = 2.*p.*(1-p);
pTwoSpots = p.^2;

%% Count the observed number of nuclei with 0,1 and 2 spots in each frame

% initialize arrays to store data
DataZeroSpots = [];
DataOneSpot = [];
DataTwoSpots = [];
AreaThreshold = 10000;

for frame = 1:length(NucleiPerFrame)
    
    Nuclei = NucleiPerFrame(frame).Nuclei;
    NThisFrame = length(Nuclei);
    
    ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold);
    
    [FractionZeroSpots, FractionOneSpot, ~, FractionTwoSpots,~, ~,TotalNuclei] = ...
        NumberOfSpotsFractions(ParticlesPerNucleus,fractionCompetent);
    
   % store the counts
    DataZeroSpots(frame) = round(FractionZeroSpots*NThisFrame);
    DataOneSpot(frame) = round(FractionOneSpot*NThisFrame);
    DataTwoSpots(frame) = round(FractionTwoSpots*NThisFrame);
end

%% Bootstrap the nmber of nuclei with 0, 1 or 2 spots
%cObsP = @(x) (sum(x)/length(x)); %bootstrapped function: fraction
cObsP = @(x) sum(x); %bootstrapped function: counts
nSamples = 1000;

ErrorZero = [];
ErrorOne = [];
ErrorTwo = [];

for frame = 1:length(NucleiPerFrame)
    %this is are the original samples in a 'bootsrappable' format
    NZeroSpot = [ones(1,DataZeroSpots(frame)) zeros(1,N(frame)-DataZeroSpots(frame))]; 
    NOneSpot = [ones(1,DataOneSpot(frame)) zeros(1,N(frame)-DataOneSpot(frame))]; 
    NTwoSpot = [ones(1,DataTwoSpots(frame)) zeros(1,N(frame)-DataTwoSpots(frame))]; 
    
    %generate the bootsrapped distributions
    [bootObsNZero,~] = bootstrp(nSamples, cObsP,NZeroSpot); 
    [bootObsNOne,~] = bootstrp(nSamples, cObsP,NOneSpot); 
    [bootObsNTwo,~] = bootstrp(nSamples, cObsP,NTwoSpot); 
    
    %calculate the bootstrapped standard error of the mean
    ErrorZero(frame) = std(bootObsNZero);
    ErrorOne(frame) = std(bootObsNOne);
    ErrorTwo(frame) = std(bootObsNTwo);
end



%% Compare prediction with observation
figure
hold on
plot(pZeroSpots .* N,'k-o','LineWidth',3)
plot(DataZeroSpots,'k-','LineWidth',3)

plot(pOneSpot .* N,'b-o','LineWidth',3)
plot(DataOneSpot,'b-','LineWidth',3)

plot(pTwoSpots .* N,'g-o','LineWidth',3)
plot(DataTwoSpots,'g-','LineWidth',3)

hold off

%% compare prediction with observation as bar graphs

ExpectedZero = (pZeroSpots.*N);
PointExpectedZero = mean(ExpectedZero(pickedFrame-1:pickedFrame+1));
ExpectedOne = (pOneSpot.*N);
PointExpectedOne = mean(ExpectedOne(pickedFrame-1:pickedFrame+1));
ExpectedTwo = (pTwoSpots.*N);
PointExpectedTwo = mean(ExpectedTwo(pickedFrame-1:pickedFrame+1));

ObservedZero = mean(DataZeroSpots((pickedFrame-1:pickedFrame+1)));
ObservedOne = mean(DataOneSpot((pickedFrame-1:pickedFrame+1)));
ObservedTwo = mean(DataTwoSpots((pickedFrame-1:pickedFrame+1)));

PointErrorZero = mean(ErrorZero((pickedFrame-1:pickedFrame+1)));
PointErrorOne = mean(ErrorOne((pickedFrame-1:pickedFrame+1)));
PointErrorTwo = mean(ErrorTwo((pickedFrame-1:pickedFrame+1)));


figure
DataForBars = [PointExpectedZero ObservedZero;PointExpectedOne ObservedOne;PointExpectedTwo ObservedTwo];
bar(DataForBars)
legend({'expectation','observation'})


figure
hold on
errorbar(1,ObservedZero,PointErrorZero,'ko','LineWidth',2,'CapSize',0)
plot(1,PointExpectedZero,'ro','MarkerFaceColor','r','MarkerSize',10)

errorbar(2,ObservedOne,PointErrorOne,'ko','LineWidth',2,'CapSize',0)
plot(2,PointExpectedOne,'ro','MarkerFaceColor','r','MarkerSize',10)

errorbar(3,ObservedTwo,PointErrorTwo,'ko','LineWidth',2,'CapSize',0)
plot(3,PointExpectedTwo,'ro','MarkerFaceColor','r','MarkerSize',10)

hold off
xlim([0 4])
ylim([0 max([ObservedZero,ObservedOne,ObservedTwo,PointExpectedZero,PointExpectedOne,PointExpectedTwo])*1.5])


%% compare prediction vs observation as time traces
AbsTime = [FrameInfo.Time]./60;
figure
subplot(1,3,1)
errorbar(AbsTime,DataZeroSpots,ErrorZero,'mo','CapSize',0,'MarkerFaceColor','m',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold on
plot(AbsTime,ExpectedZero,'o','MarkerFaceColor','k',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold off

subplot(1,3,2)
errorbar(AbsTime,DataOneSpot,ErrorOne,'go','CapSize',0,'MarkerFaceColor','g',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold on
plot(AbsTime,ExpectedOne,'o','MarkerFaceColor','y',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold off

subplot(1,3,3)
errorbar(AbsTime,DataTwoSpots,ErrorTwo,'bo','CapSize',0,'MarkerFaceColor','b',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold on
plot(AbsTime,ExpectedTwo,'o','MarkerFaceColor','c',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold off

%% just show the observed number of nuclei with zero, one or two spots

figure
hold on
errorbar(AbsTime,DataZeroSpots,ErrorZero,'mo','CapSize',0,'MarkerFaceColor','m',...
    'MarkerEdgeColor','none','MarkerSize',10)
errorbar(AbsTime,DataOneSpot,ErrorOne,'go','CapSize',0,'MarkerFaceColor','g',...
    'MarkerEdgeColor','none','MarkerSize',10)
errorbar(AbsTime,DataTwoSpots,ErrorTwo,'bo','CapSize',0,'MarkerFaceColor','b',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold off
legend({'zero','one','two'})
xlabel('time (min)')

%% Bootstrap the FRACTION of nuclei with 0, 1 or 2 spots
cObsP = @(x) (sum(x)/length(x)); %bootstrapped function: fraction
%cObsP = @(x) sum(x); %bootstrapped function: counts
nSamples = 1000;

ErrorZero = [];
ErrorOne = [];
ErrorTwo = [];

for frame = 1:length(NucleiPerFrame)
    %this is are the original samples in a 'bootsrappable' format
    NZeroSpot = [ones(1,DataZeroSpots(frame)) zeros(1,N(frame)-DataZeroSpots(frame))]; 
    NOneSpot = [ones(1,DataOneSpot(frame)) zeros(1,N(frame)-DataOneSpot(frame))]; 
    NTwoSpot = [ones(1,DataTwoSpots(frame)) zeros(1,N(frame)-DataTwoSpots(frame))]; 
    
    %generate the bootsrapped distributions
    [bootObsFZero,~] = bootstrp(nSamples, cObsP,NZeroSpot); 
    [bootObsFOne,~] = bootstrp(nSamples, cObsP,NOneSpot); 
    [bootObsFTwo,~] = bootstrp(nSamples, cObsP,NTwoSpot); 
    
    %calculate the bootstrapped standard error of the mean
    ErrorFZero(frame) = std(bootObsFZero);
    ErrorFOne(frame) = std(bootObsFOne);
    ErrorFTwo(frame) = std(bootObsFTwo);
end


%% show the FRACTION of nuclei with zero, one or two spots
%close all
figure
hold on
errorbar(AbsTime,DataZeroSpots./N,ErrorFZero,'mo','CapSize',0,'MarkerFaceColor','m',...
    'MarkerEdgeColor','none','MarkerSize',10)
errorbar(AbsTime,DataOneSpot./N,ErrorFOne,'go','CapSize',0,'MarkerFaceColor','g',...
    'MarkerEdgeColor','none','MarkerSize',10)
errorbar(AbsTime,DataTwoSpots./N,ErrorFTwo,'bo','CapSize',0,'MarkerFaceColor','b',...
    'MarkerEdgeColor','none','MarkerSize',10)
hold off
legend({'zero','one','two'})
xlabel('time (min)')
ylabel('fraction of nuclei')
ylim([0 1])




end
