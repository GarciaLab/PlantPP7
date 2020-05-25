function [bootstrappedMean,bootstrappedStd] = bootstrapFracNSpots(Ntot,Ni)

% returns the bootstrapped mean and standard deviation in the fraction of
% actively transcribing nuclei
% Arguments: 
% Ntot is the total number of nuclei in the sample
% Ni is the observed number of nuclei with i spots (i = {0,1,2})

nSamples = 10000; %how many times we're bootstrapping
cObsP = @(x) (sum(x)/length(x)); %this is the function we are bootstrapping, a ratio

% bootstrapping the observed fraction of one spot nuclei
NoneSpot = [ones(1,Ni) zeros(1,Ntot-Ni)]; %this is the original sample
[bootObsNOne,~] = bootstrp(nSamples, cObsP,NoneSpot); 

bootstrappedMean = mean(bootObsNOne);
bootstrappedStd = std(bootObsNOne);

end


