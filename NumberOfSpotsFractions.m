function [FractionZeroSpots, FractionOneSpot, FractionOneOrTwoSpots, FractionTwoSpots,...
    FractionTwoOrMoreSpots, FractionActiveNuclei,TotalNuclei] = ...
    NumberOfSpotsFractions (ParticlesPerNucleus)

[N,~] = histcounts(ParticlesPerNucleus);
% sum(N) is total nuclei, N(1) is nuclei with 0 spots, N(2) is nuclei with one, N(3) have two and
% sum(N(3:end)) is two or more.

FractionZeroSpots = N(1)/sum(N);

if length(N) > 1
    FractionOneSpot = N(2)/sum(N);
else
    FractionOneSpot =0;
end

if length(N)>2
    FractionOneOrTwoSpots = (N(2)+N(3))/sum(N);
    FractionTwoSpots = N(3)/sum(N);
else
    FractionOneOrTwoSpots = 0;
    FractionTwoSpots = 0;
end

FractionTwoOrMoreSpots = [];%(sum(N)-N(1)-N(2))/sum(N);
FractionActiveNuclei = (sum(N)-N(1))/sum(N);
TotalNuclei = sum(N);

end
