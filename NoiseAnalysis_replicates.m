%%
FirstSpot = [load('9.4rep2_spots1_integrated.mat')...
    load('9.4rep1_spots1_integrated.mat')];
FirstSpots = [FirstSpot.RawSpots1];

SecondSpot = [load('9.4rep2_spots2_integrated.mat')...
    load('9.4rep1_spots2_integrated.mat')];
SecondSpots = [SecondSpot.RawSpots2];

% scramble the data for visualization purposes
[FirstSpots,SecondSpots] = scrambleRows(FirstSpots,SecondSpots);

MeanAll = mean([FirstSpots SecondSpots]);
NormFirstSpots = FirstSpots./MeanAll;
NormSecondSpots = SecondSpots./MeanAll;

Spots1Fluos = 1-NormFirstSpots;
Spots2Fluos = 1-NormSecondSpots;
figure
plot(NormFirstSpots,NormSecondSpots,'ko')
hold on
plot([0 max([NormFirstSpots NormSecondSpots])*1.2],...
    [0 max([NormFirstSpots NormSecondSpots])*1.2],'k-')
title('normalized integrated fluorescence')
% bootstrap the errors to get errobars
figure
subplot(1,2,1)
plot(FirstSpots,SecondSpots,'o')
hold on
plot([0 max([FirstSpots SecondSpots])*1.2],...
    [0 max([FirstSpots SecondSpots])*1.2],'k-')

nSamples = 10000; %how many times we're bootstrapping
%these are the function we are bootstrapping, the formulas for each noise component
TotObsP = @(x,y) 1/2 * (mean(x.^2) + mean(y.^2)); 
IntObsP = @(x,y) 1/2 * mean((x-y).^2);
CorrObsP = @(x,y) mean(x.*y);

% bootstrapping the observed total noise
[bootTotObs,~] = bootstrp(nSamples,TotObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanTotalNoise = mean(bootTotObs);
bootstrappedStdTotalNoise = std(bootTotObs);

% bootstrapping the observed intrinsic noise
[bootIntObs,~] = bootstrp(nSamples,IntObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanIntrinsicNoise = mean(bootIntObs);
bootstrappedStdIntrinsicNoise = std(bootIntObs);

% bootstrapping the observed correlated noise
[bootCorrObs,~] = bootstrp(nSamples,CorrObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanCorrNoise = mean(bootCorrObs);
bootstrappedStdCorrNoise = std(bootCorrObs);

subplot(1,2,2)
hold on
errorbar(1,bootstrappedMeanTotalNoise,bootstrappedStdTotalNoise,'bo','LineWidth',3,...
    'CapSize',0)
errorbar(2,bootstrappedMeanIntrinsicNoise,bootstrappedStdIntrinsicNoise,'ro','LineWidth',3,...
    'CapSize',0)
errorbar(3,bootstrappedMeanCorrNoise,bootstrappedStdCorrNoise,'go','LineWidth',3,...
    'CapSize',0)
hold off
ylim([0 1])
xlim([0.8 3.2])
legend('\eta_{tot}^2','\eta_{int}^2','\eta_{corr}^2')