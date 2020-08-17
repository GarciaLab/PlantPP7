function TwoSpotSnapshots_v2(Prefix,FractionCompetent)
%takes a single Prefix and generates figures related to two spots

%% load nuclear mask, particles.mat and spots.mat
% this is a single binary image where 1s correspond to nuclei and 0s are
% background
DynamicsResultsPath = '/Users/simon_alamos/Dropbox/DynamicsResults';
PrefixPath = [DynamicsResultsPath '/' Prefix];
NuclearMaskPath = [PrefixPath '/FilledFilteredMask.mat'];
try
    load(NuclearMaskPath);
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

%% flag overlapping particles so that we don't count them twice
Particles = disapproveOverlappingParticles(Particles,Spots);

%% Assign particles to nuclei from scratch
% this script uses CompiledParticles, make sure to recompile after
% disapproving overlapping ones
% CompiledParticles = CompileParticles(Prefix);
Nuclei = assignParticlesToNuclei (FilledFilteredMask,CompiledParticles);

% Particles per nucleus: analysis of results
AreaThreshold = 6000; % maximum number of pixels squared for a nucleus to be considered
% nuclei larger than this are polyploid
% now make a vector containing the number of spots per nucleus
ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold);

%% Count the number of nuclei with zero, one or two spots

[FractionZeroSpots, FractionOneSpot, FractionOneOrTwoSpots, FractionTwoSpots,...
    FractionTwoOrMoreSpots, FractionActiveNuclei,TotalNuclei] = NumberOfSpotsFractions (ParticlesPerNucleus,1);

FractionOneSpot = FractionOneSpot / FractionCompetent;
FractionTwoSpots = FractionTwoSpots / FractionCompetent;

%% Calculate the bootstrapped number of observed two spot nuclei and the independence prediction

[bootstrappedOneSpotError,bootstrappedTwoSpotError,...
ExpectedTwoSpot,ExpectedErrorTwoSpot,...
Predicted_p,ErrorPredicted_p] = ...
MultinomialTestOfTwoSpots_time_v2(FractionOneSpot,FractionTwoSpots,ceil(TotalNuclei*FractionCompetent));


% Plot measurement vs prediction side by side

% now plot this as bar graphs for a single frame
% the prediction
Prediction = Predicted_p.^2;

figure
errorbar(0.9,Prediction,abs(Prediction-ExpectedErrorTwoSpot(1)),abs(Prediction-ExpectedErrorTwoSpot(2)),'ko')
hold on
% the bootstrapped measurement
errorbar(1,FractionTwoSpots,bootstrappedTwoSpotError,'bo')
errorbar(1.1,FractionOneSpot,bootstrappedOneSpotError,'go')
hold off
xlim([0.8 1.2])
ylim([0 0.4])
legend('prediction','measurement')
title(['Measured vs Predicted Two Spots, fraction competent = ' num2str(FractionCompetent) '_' Prefix],'interpreter','none')
savefig([PrefixPath '/TwoSpots_bars_fraction_competent.fig'])






end
