function MeanInstaFractionOn = getMeanInstaFractionOn(Struct)

%find the longest movie to define the last time point
Durations = [];
for rep = 1:length(Struct)
    repAbsTime = Struct(rep).AbsTime;
    Durations = [Durations repAbsTime(end)];
end
LastTimePoint = ceil(max(Durations));

fractionOns = nan(rep,LastTimePoint);


for t = 1:LastTimePoint
    for rep = 1:length(Struct)
        repAbsTime = Struct(rep).AbsTime;
        % get the right particles of this rep in this time point
        [dummy,nearestFrame] = min(abs(Struct(rep).AbsTime-t)); 
        thisRepFracOnThisTime = Struct(rep).InstFractionON;
        thisRepFracOnThisTime = thisRepFracOnThisTime(nearestFrame);
        
        fractionOns(rep,t) = thisRepFracOnThisTime;
    end
end

%errorbar(nanmean(fractionOns),std(fractionOns)./sqrt(rep))

MeanInstaFractionOn = nanmean(fractionOns);



end

