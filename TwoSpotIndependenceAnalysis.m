function TwoSpotIndependenceAnalysis(Prefix,DynamicsResultsPath,pickedFrame,fractionCompetent)


%% Load stuff
PrefixPath = [DynamicsResultsPath '/' Prefix];
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

% initialize arrays to store data
DataOneSpot =[];
DataTwoSpots = [];
ExpectedOneSpot = [];
ExpectedTwoSpot = [];
ExpectedErrorOneSpot = [];
ExpectedErrorTwoSpot = nan(2,length(NucleiPerFrame));
Predicted_p = [];

%%
% In each frame we will use the fraction of nuclei with zero and one spot to predict the
% fraction of nuclei with two spots. Then we'll compare it to data.
FractionCompetent = 1;
TotalNucleiPerFrame = [];
AreaThreshold = 10; %pixels
for frame = 1:length(NucleiPerFrame)
    % count stuff
    Nuclei = NucleiPerFrame(frame).Nuclei;
    ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold);
    [~, FractionOneSpot, ~, FractionTwoSpots,~, ~,TotalNuclei] = ...
        NumberOfSpotsFractions(ParticlesPerNucleus,fractionCompetent);
    
    % store the counts
    DataOneSpot(frame) = FractionOneSpot*TotalNuclei;
    DataTwoSpots(frame) = FractionTwoSpots*TotalNuclei;
    TotalNucleiPerFrame = [TotalNucleiPerFrame TotalNuclei];
    
    [bootstrappedOneSpotError(frame),bootstrappedTwoSpotError(frame),...
        ExpectedTwoSpot(frame),ExpectedErrorTwoSpot(:,frame),...
        Predicted_p(frame),ErrorPredicted_p(frame)] = ...
    MultinomialTestOfTwoSpots_time_v2(FractionOneSpot,FractionTwoSpots,TotalNuclei*FractionCompetent);
end

%% plots

figure
hold on
%plot fitted p
shadedErrorBar(AbsTime,Predicted_p,ErrorPredicted_p)
%plot data
shadedErrorBar(AbsTime,DataOneSpot./(TotalNucleiPerFrame.*FractionCompetent),...
    bootstrappedOneSpotError,'lineProps',{'Color','g','LineWidth',2})
ylabel('fraction of nuclei and  p')
xlabel('time (min)')
hold off
title(['One Spot ' Prefix],'interpreter','none')
savefig([PrefixPath '/OneSpot_time.fig'])


figure
hold on
%plot predicted fraction with two spots
shadedErrorBar(AbsTime,Predicted_p.^2,...
    flipud(ExpectedErrorTwoSpot))%./(TotalNucleiPerFrame.*FractionCompetent))
%plot measured observations
shadedErrorBar(AbsTime,DataTwoSpots./(TotalNucleiPerFrame.*FractionCompetent),...
    bootstrappedTwoSpotError,'lineProps',{'Color','b','LineWidth',2})
ylabel('fraction of nuclei with two spots and p^2')
xlabel('time (min)')
title(['Two Spot ' Prefix],'interpreter','none')
hold off
savefig([PrefixPath '/TwoSpots_time.fig'])


% now plot this as bar graphs for a single frame
ConfidenceBounds = flipud(ExpectedErrorTwoSpot)./(TotalNucleiPerFrame.*FractionCompetent);
% the prediction
Prediction = Predicted_p.^2;

figure
errorbar(0.9,Prediction(pickedFrame),ConfidenceBounds(2,pickedFrame),ConfidenceBounds(1,pickedFrame),'o')
hold on
% the bootstrapped measurement
Data = DataTwoSpots./(TotalNucleiPerFrame.*FractionCompetent);
errorbar(1.1,Data(pickedFrame),bootstrappedTwoSpotError(pickedFrame),'o')
hold off
xlim([0.8 1.2])
legend('prediction','measurement')
title(['Measured vs Predicted Two Spots ' Prefix],'interpreter','none')
savefig([PrefixPath '/TwoSpots_bars.fig'])








% 
% % First, the prediction:
% % bootstrap the number of nuclei with one spot to claculate a spread in p^2
% NNuclei = round(mean(TotalNucleiPerFrame(firstFrame:lastFrame)));
% MeasuredOneS = round(mean(DataOneSpot(firstFrame:lastFrame)));
% cObsP = @(x) (sum(x)/length(x)); %bootstrapped function: fraction
% nSamples = 1000;
% NoOneS = [ones(1,MeasuredOneS) zeros(1,NNuclei-MeasuredOneS)]; %this is the original sample
% [bootObsNOneS,~] = bootstrp(nSamples, cObsP,NoOneS); % this is the bootstrapped sample
% bootObsNOneS(bootObsNOneS>0.5) = nan; %
% bootstrapped_p = 0.5*(1 - sqrt(1-(2.*bootObsNOneS))); %bootstrapped p^2
% bootstrapped_p_sqrd = bootstrapped_p.^2;
% PredictedTwoS = nanmean(bootstrapped_p_sqrd);
% 
% % Here we deal with the error in the prediction
% Var_two_spots = NNuclei.*bootstrapped_p_sqrd.*(1-bootstrapped_p_sqrd);
% SD_two_spots = sqrt(Var_two_spots);
% % calculate the confidence interval without assuming normality
% CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
% ExpectedErrorTwoSpot = CIFcn(SD_two_spots,95)./NNuclei; %returns lower and upper interval values in the number of nuclei
% 
% % Now, the measurement. 
% % bootstrap the number of nuclei with two spots
% MeasuredTwoS = round(mean(DataTwoSpots(firstFrame:lastFrame)));
% NoTwoS = [ones(1,MeasuredTwoS) zeros(1,NNuclei-MeasuredTwoS)]; %this is the original sample
% [bootObsNOneS,~] = bootstrp(nSamples, cObsP,NoTwoS);
% errorMeasuredTwoS = std(bootObsNOneS);
% 
% 
% figure
% errorbar([0.1 0.2],[MeasuredOneS/NNuclei PredictedTwoS],...
%     [errorMeasuredTwoS ExpectedErrorTwoSpot(1)],[errorMeasuredTwoS ExpectedErrorTwoSpot(2)],'o')
% xlim([0 0.3])
% title(['Measured vs Predicted Two Spots ' Prefix],'interpreter','none')
% savefig([PrefixPath '/TwoSpots_bars.fig'])


end
