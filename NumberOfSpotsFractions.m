function [FractionZeroSpots, FractionOneSpot, FractionOneOrTwoSpots, FractionTwoSpots,...
    FractionTwoOrMoreSpots, FractionActiveNuclei] = ...
    NumberOfSpotsFractions (ParticlesPerNucleus)

[N,~] = histcounts(ParticlesPerNucleus);
% sum(N) is total nuclei, N(1) is nuclei with 0 spots, N(2) is nuclei with one, N(3) have two and
% sum(N(3:end)) is two or more.

FractionZeroSpots = N(1)/sum(N);
FractionOneSpot = N(2)/sum(N);
FractionOneOrTwoSpots = (N(2)+N(3))/sum(N);
FractionTwoSpots = N(3)/sum(N);
FractionTwoOrMoreSpots = (sum(N)-N(1)-N(2))/sum(N);
FractionActiveNuclei = (sum(N)-N(1))/sum(N);

end
