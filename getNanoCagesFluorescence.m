function [nanocagesFluorescenceVector] = getNanoCagesFluorescence(Prefixes,PathToDynamicsResults)


MaxSpotNumber = 1000; %this number doesn't matter, it just has to be larger than the number of datapoints (i.e nanocages in all the datasets)

% AllFixedAreaIntensity_60 = nan(length(Prefixes),MaxSpotNumber);
% Allgauss3DIntensity_60 = nan(length(Prefixes),MaxSpotNumber);
% Allgauss3DIntensityVector_60 = [];
nanocagesFluorescenceVector =[];

for p = 1:length(Prefixes)
    Prefix = Prefixes{p};
    PrefixPath = [PathToDynamicsResults '/' Prefix];
    load([PrefixPath '/Spots.mat']);
    AllSpots = Spots.Fits;
    MaxFixedAreaIntensity = [];

    for s = 1:length(AllSpots)
        MaxFixedAreaIntensity(s) = max(AllSpots(s).FixedAreaIntensity);
    end
    
    nanocagesFluorescenceVector = [nanocagesFluorescenceVector MaxFixedAreaIntensity];
end
