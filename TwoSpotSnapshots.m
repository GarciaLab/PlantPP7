function TwoSpotSnapshots(Prefix)
%takes a single Prefix and generates figures related to two spots

%% load nuclear mask, particles.mat and spots.mat
% this is a single binary image where 1s correspond to nuclei and 0s are
% background
DynamicsResultsPath = '/Users/simon_alamos/Dropbox/DynamicsResults';
PrefixPath = [DynamicsResultsPath '/' Prefix];

NuclearMaskPath = [PrefixPath '/FilledFilteredMask.mat'];
try
    load(NuclearMaskPath)
catch
end
    
load([PrefixPath '/' 'CompiledParticles.mat'],'CompiledParticles')
%in some versions of the code CompiledParticles can be a struct within a cell
if iscell(CompiledParticles) 
    CompiledParticles = CompiledParticles{1};
end

load([PrefixPath '/' 'Spots.mat'])
load([PrefixPath '/' 'Particles.mat'],'Particles')
load([PrefixPath '/' 'FrameInfo.mat'])
load([PrefixPath '/' 'FrameInfo.mat'])

try
load([PrefixPath '/' 'Ellipses.mat'])
catch
end
%% In the case of time lapse movies run this part only

[NucleiPerFrame, CompiledParticles] = ...
    assignParticlesToEllipses(Ellipses,CompiledParticles,FrameInfo);

figure
hold on
for frame = 1:length(NucleiPerFrame)
    Nuclei = NucleiPerFrame(frame).Nuclei;
    ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,10);
    [~, FractionOneSpot, ~, FractionTwoSpots,~, ~,TotalNuclei] = ...
        NumberOfSpotsFractions (ParticlesPerNucleus);
    MultinomialTestOfTwoSpots_norm(FractionOneSpot,FractionTwoSpots,ceil(TotalNuclei*0.5),0)
end
hold off
title(Prefix,'interpreter','none')
savefig([PrefixPath '/OneTwoSpots_scatter.fig'])
ylim([0 1])

DataOneSpot =[];
DataTwoSpots = [];
ExpectedOneSpot = [];
ExpectedTwoSpot = [];
ExpectedErrorOneSpot = [];
ExpectedErrorTwoSpot = [];

for frame = 1:length(NucleiPerFrame)
    
    Nuclei = NucleiPerFrame(frame).Nuclei;
    ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,10);
    [~, FractionOneSpot, ~, FractionTwoSpots,~, ~,TotalNuclei] = ...
        NumberOfSpotsFractions (ParticlesPerNucleus);
    
    DataOneSpot(frame) = FractionOneSpot*TotalNuclei;
    DataTwoSpots(frame) = FractionTwoSpots*TotalNuclei;
    
    [ExpectedOneSpot(frame),ExpectedTwoSpot(frame),...
        ExpectedErrorOneSpot(frame),ExpectedErrorTwoSpot(frame)] = ...
    MultinomialTestOfTwoSpots_time(FractionOneSpot,FractionTwoSpots,TotalNuclei,frame)
end


AbsTime = [FrameInfo.Time]./60; %in minutes
figure
hold on
shadedErrorBar(AbsTime,ExpectedOneSpot,ExpectedErrorOneSpot,'lineProps',...
    {'Color','b','LineWidth',2,'LineStyle','--'})
plot(AbsTime,DataOneSpot,'b-','LineWidth',2,'MarkerFaceColor','b')
title(['One Spot ' Prefix],'interpreter','none')
hold off
savefig([PrefixPath '/OneSpot_time.fig'])
ylim([0 max(ExpectedOneSpot+ExpectedErrorOneSpot)*1.2])


figure
hold on
shadedErrorBar(AbsTime,ExpectedTwoSpot,ExpectedErrorTwoSpot,'lineProps',...
    {'Color','r','LineWidth',2,'LineStyle','--'})
plot(AbsTime,DataTwoSpots,'r-','LineWidth',2,'MarkerFaceColor','r')
title(['Two Spot ' Prefix],'interpreter','none')
hold off
savefig([PrefixPath '/TwoSpots_time.fig'])
ylim([0 max(DataTwoSpots)*1.2])

%% flag overlapping particles so that we don't count them twice
Particles = disapproveOverlappingParticles(Particles,Spots);

%% Assign particles to nuclei from scratch
% this script uses CompiledParticles, make sure to recompile after
% disapproving overlapping ones
% CompiledParticles = CompileParticles(Prefix);
Nuclei = assignParticlesToNuclei (FilledFilteredMask,CompiledParticles);

%% Particles per nucleus: analysis of results
AreaThreshold = 6000; % maximum number of pixels squared for a nucleus to be considered
% nuclei larger than this are polyploid
% now make a vector containing the number of spots per nucleus
ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold);

[FractionZeroSpots, FractionOneSpot, FractionOneOrTwoSpots, FractionTwoSpots,...
    FractionTwoOrMoreSpots, FractionActiveNuclei,TotalNuclei] = NumberOfSpotsFractions (ParticlesPerNucleus);
close all
hold on
MultinomialTestOfTwoSpots_norm(FractionOneSpot,FractionTwoSpots,63,1)
hold off
title(Prefix,'interpreter','none')
savefig([PrefixPath '/OneTwoSpots_scatter.fig'])



%%


% P = probability of single locus turning on

% calculate P based on the fraction off, under assumption of independence
% the expected probability that both are off is (1-P)^2, solving for P
% using the measured  fraction of cells with both loci off:
PfromPoff = 1-sqrt(FractionZeroSpots); 

% calculate P based on the fraction with one or two, under assumption of independence
% the expected probability that both are off is (1-P)^2, solving for P
% using the measured  fraction of cells with both loci off:
PfromOneOrTwo = 1-sqrt(1-FractionOneOrTwoSpots);

% now predict the fraction with two spots, will take the average of the two
% predictions
PredictedTwoSpots = (mean([PfromPoff PfromOneOrTwo]))^2;
% prediction of one spot
PredictedOneSpot = 2*mean([PfromPoff PfromOneOrTwo])*(1-mean([PfromPoff PfromOneOrTwo]));

nSamples = 10000;
% bootstrapping the onserved fraction of off nuclei
NzeroSpots = [ones(1,OffNuclei) zeros(1,TotalNuclei-OffNuclei)]; %this is the original sample
cObsP = @(x) (sum(x)/length(x));
[bootObsNoff,~] = bootstrp(nSamples, cObsP, NzeroSpots);
% mean(bootstat) is the bootstrapped mean number of off nuclei
% std(bootstat) is the standard error of the mean number of off nuclei
figure
histogram(bootObsNoff)
title('bootstrapped observed zero spots')

% bootstrapping the observed fraction of one spot nuclei
NoneSpot = [ones(1,OneSpot) zeros(1,TotalNuclei-OneSpot)]; %this is the original sample
[bootObsNOne,~] = bootstrp(nSamples, cObsP,NoneSpot); 
% mean(bootstat) is the bootstrapped mean number of off nuclei
% std(bootstat) is the standard error of the mean number of off nuclei
figure
histogram(bootObsNOne)
title('bootstrapped observed one spot')

% bootstrapping the observed number two spot nuclei
NtwoSpot = [ones(1,TwoOrMore) zeros(1,TotalNuclei-TwoOrMore)]; %this is the original sample
[bootObsNTwo,~] = bootstrp(nSamples, cObsP,NtwoSpot); 
% mean(bootstat) is the bootstrapped mean number of off nuclei
% std(bootstat) is the standard error of the mean number of off nuclei
figure
histogram(bootObsNTwo)
title('bootstrapped observed two spots')

% bootstrapping the PREDICTED probability of one spot nuclei, 2p(1-p)
c1p = @(x) 2*((1-sqrt(sum(x)/length(x))))*(1-((1-sqrt(sum(x)/length(x)))));
[bootOneSpot,~] = bootstrp(nSamples, c1p, NzeroSpots);
%bootPredOnespot = mean(bootOneSpot);
%bootPredOnespotError = std(bootOneSpot);% is the standard error of the calculated p
figure
histogram(bootOneSpot)
title('bootstrapped predicted one spot')

% bootstrapping the PREDICTED probability of two spot nuclei, p^2
c2p = @(x) (1-sqrt(sum(x)/length(x))).^2;
[bootTwoSpot,~] = bootstrp(nSamples, c2p, NzeroSpots);
%bootrappedP = mean(bootTwoSpot);
%bootrappedPError = std(bootTwoSpot);% is the standard error of the calculated p
figure
histogram(bootTwoSpot)
title('bootstrapped predicted two spots')

figure % with error for each group, off, 1spot and twospot, measured and predicted
% measured values
meanOneSpotObserved = mean(bootObsNOne);
ErrorOneSpotObserved = std(bootObsNOne);
meanTwoSpotObserved= mean(bootObsNTwo);
ErrorTwoSpotObserved = std(bootObsNTwo);

% predicted values
meanOneSpotPredicted = mean(bootOneSpot);
ErrorOneSpotPredicted = std(bootOneSpot);
meanTwoSpotPredicted= mean(bootTwoSpot);
ErrorTwoSpotPredicted = std(bootTwoSpot);

hold on
errorbar(.5,meanOneSpotObserved,ErrorOneSpotObserved,'go','MarkerFaceColor','g')
errorbar(1,meanTwoSpotObserved,ErrorTwoSpotObserved,'yo','MarkerFaceColor','y')
errorbar(1.5,meanOneSpotPredicted,ErrorOneSpotPredicted,'bo','MarkerFaceColor','b')
errorbar(2,meanTwoSpotPredicted,ErrorTwoSpotPredicted,'ro','MarkerFaceColor','r')
hold off
ylim([0 1])
xlim([0 2.5])

%%
blackMask = zeros(size(NucleusLabels));
spotDensityImage = blackMask;
spotsMaskAll = blackMask;
for p = 1:length(CompiledParticles)
    particleXPos = CompiledParticles(p).xPos;
    particleYPos = CompiledParticles(p).yPos;
%     spotMask = blackMask;
%     spotMask(particleYPos,particleXPos) = 1;
    spotsMaskAll(particleYPos,particleXPos) = 1;
%     spotDensityImage = spotDensityImage + bwdist(spotMask);
     p
end

nucleiDensityImage = blackMask;
nucleiMaskAll = blackMask;
for n = 1:length(Nuclei)
    nucleusYPos = floor(Nuclei(n).Centroid(1));
    nucleusXPos = floor(Nuclei(n).Centroid(2));
%     nucleusMask = blackMask;
%     nucleusMask(nucleusYPos,nucleusXPos) = 1;
    nucleiMaskAll(nucleusYPos,nucleusXPos) = 1;
%     nucleiDensityImage = nucleiDensityImage + bwdist(nucleusMask);
    n
end

%imshow(log10(spotDensityImage./nucleiDensityImage),[])
% colormap viridis
% N = bwdist(nucleiMaskAll);
% P = bwdist(spotsMaskAll);
% D = P./imcomplement(N);
% imshow(D,[])


windowWidth = 5; % Whatever you want.  More blur for larger numbers.
%kernel = ones(windowWidth) / windowWidth ^ 2;
H = fspecial('gaussian',2000,300);
%blurredSpotsImage = imfilter(spotsMaskAll, kernel); % Blur the image.
blurredSpotsImage = imfilter(spotsMaskAll, H)+1;
% blurredNucleiImage = imfilter(nucleiMaskAll, kernel);
blurredNucleiImage = imfilter(nucleiMaskAll, H)+1;
D = blurredSpotsImage./blurredNucleiImage;
imshow(blurredSpotsImage./blurredNucleiImage,[]); % Display it.
colormap viridis

%% Voronoi cells
XY = ceil([Nuclei.Centroid]);
Y = XY(1:2:end);
X = XY(2:2:end);
[vx,vy] = voronoi(X,Y);

[v c] = voronoin([X',Y'],{'Qbb'});

figure
Palette =[parula(8)+white(8)]./2;
hold on
for i = 1:length(c)
    NumberOfSpots = Nuclei(i).NParticles;
    if NumberOfSpots == 0
        Color = Palette(1,:);
    elseif NumberOfSpots == 1
        Color = Palette(5,:);
    else
        Color = Palette(end,:);
    end
    fill(v(c{i},1),v(c{i},2),Color,'LineWidth',1,'EdgeColor','k')
end
hold off
ylim([0 4096])
xlim([0 4096])
set(gca,'color',[0 0 0])

%% number of particles as a function of nucleus area
figure
scatter([Nuclei.Area],ParticlesPerNucleus,'ro')
xlabel('nucleus area (pixels)^3')
ylabel('number of spots per nucleus')
%% histogram of nucleus area for each number of particles
figure
NumberOfSpots = [Nuclei.NumberOfSpots]
NucleiAreas = sqrt([Nuclei.Area]).^3;
hold on
for NSpots = 0:max(NumberOfSpots)
    BinNucleiIDs = NumberOfSpots == NSpots;
    BinNucleiAreas = NucleiAreas(BinNucleiIDs);
    histogram(BinNucleiAreas,10)
end
hold off

%% Distribution of particle fluorescence
FluoVals = [CompiledParticles.Fluo];
FluoVals = FluoVals(FluoVals>0);
histogram(log(FluoVals),30)

%% spatial distribution on/off nuclei
ColorMap = [hsv(5) + white(5)]./2;

OnOffNucleiLabels = DilatedNucleusLabels;
for n = 1:length(Nuclei)
    if Nuclei(n).NParticles == 1    
        OnOffNucleiLabels(OnOffNucleiLabels==n)=2;
    elseif Nuclei(n).NParticles == 2
        OnOffNucleiLabels(OnOffNucleiLabels==n)=5;
    else
        OnOffNucleiLabels(OnOffNucleiLabels==n)=1;
    end
end
imshow(OnOffNucleiLabels,[])
colormap

title('on and off nuclei')

%% closest neighbors fraction active
AllNucleiCentroids = [Nuclei.Centroid];
AllNucleiXPos = AllNucleiCentroids(1:2:end);
AllNucleiYPos = AllNucleiCentroids(2:2:end);
OnOffNuclei = ([Nuclei.NParticles])>0;
FractionOn = sum(OnOffNuclei)/length(OnOffNuclei);
RadialFractionActive = nan(sum(OnOffNuclei),length(Nuclei));
hold on
counter = 0;
for n = 1:length(Nuclei)
    FractionActiveNeighborhood = nan(1,length(Nuclei));
    if ~isempty(Nuclei(n).AssociatedParticles2)
        counter = counter+1;
        for k = 2:length(Nuclei)
            OnNucleusXPos = Nuclei(n).Centroid(1);
            OnNucleusYPos = Nuclei(n).Centroid(2);
            DistToRest = sqrt((OnNucleusXPos-AllNucleiXPos).^2 + (OnNucleusYPos-AllNucleiYPos).^2);
            DistToRest(DistToRest==0)=nan; %remove distance to self
            SortedDistances = sort(DistToRest);
            OnNeighbors = sum(OnOffNuclei(ismember(DistToRest,SortedDistances(1:k))));
            FractionActiveNeighborhood(k) = (OnNeighbors+1)/(k+1);
        end
        Nuclei(n).FractionActiveClosestK = FractionActiveNeighborhood./FractionOn;
        plot(Nuclei(n).FractionActiveClosestK)
        RadialFractionActive(counter,:) = FractionActiveNeighborhood./FractionOn;
    end
end
hold off
figure
errorbar(nanmean(RadialFractionActive),nanstd(RadialFractionActive,1)./sqrt(counter),'k-','CapSize',0)



%% spatial distribution spot intensities
% if a nucleus has more than one, add their intensities
SpotIntensityNucleiLabels = DilatedNucleusLabels;
for n = 1:length(Nuclei)
    if Nuclei(n).AssociatedParticles1 
        ThisNucleusParticlesFluo = nan(1,length(Nuclei(n).AssociatedParticles1));
        for p = 1:length(ThisNucleusParticlesFluo)
            particleID = Nuclei(n).AssociatedParticles1(p);
            ThisNucleusParticlesFluo(p) = CompiledParticles(particleID).Fluo;
        end
        SpotIntensityNucleiLabels(SpotIntensityNucleiLabels==n)= nansum(ThisNucleusParticlesFluo);
    else
        SpotIntensityNucleiLabels(SpotIntensityNucleiLabels==n)=1;
    end
end
figure
imshow(SpotIntensityNucleiLabels,[])
colormap plasma
title('sum of spot fluorescence per nucleus')
figure
histogram(unique(SpotIntensityNucleiLabels(SpotIntensityNucleiLabels>1)),30)
title('sum of spot fluorescence per nucleus')
xlabel('total spot intensity per nucleus (a.u)')
ylabel('counts')

%% spatial distribution spot intensities
% if a nucleus has more than one, average their intensities
SpotIntensityNucleiLabels = DilatedNucleusLabels;
for n = 1:length(Nuclei)
    if Nuclei(n).AssociatedParticles1 
        ThisNucleusParticlesFluo = nan(1,length(Nuclei(n).AssociatedParticles1));
        for p = 1:length(Nuclei(n).AssociatedParticles1)
            particleID = Nuclei(n).AssociatedParticles1(p);
            ThisNucleusParticlesFluo(p) = CompiledParticles(particleID).Fluo;
        end
        SpotIntensityNucleiLabels(SpotIntensityNucleiLabels==n)= nanmean(ThisNucleusParticlesFluo);
    else
        SpotIntensityNucleiLabels(SpotIntensityNucleiLabels==n)=1;
    end
end
figure
imshow(SpotIntensityNucleiLabels,[])
colormap plasma
title('average spot fluorescence per nucleus')

figure
histogram(unique(SpotIntensityNucleiLabels(SpotIntensityNucleiLabels>1)),30)
title('average spot fluorescence per nucleus')
xlabel('mean intensity of all spots per nucleus')


%% spatial distribution number of spots/nucleus: either 0, 1 or >1
SpotsPerNucleiLabels = DilatedNucleusLabels;
for n = 1:length(Nuclei)
    if Nuclei(n).AssociatedParticles1
        if length(Nuclei(n).AssociatedParticles1) == 1
            SpotsPerNucleiLabels(SpotsPerNucleiLabels==n)= 4.5;
        elseif length(Nuclei(n).AssociatedParticles1) == 2
            SpotsPerNucleiLabels(SpotsPerNucleiLabels==n)= 6;
        else
            SpotsPerNucleiLabels(SpotsPerNucleiLabels==n)= 8;
        end
    else
        SpotsPerNucleiLabels(SpotsPerNucleiLabels==n)=2;
    end
end
imshow(SpotsPerNucleiLabels,[])
viridis(1,:) = [1 1 1];
colormap viridis


