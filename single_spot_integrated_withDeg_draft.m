gamma = 0.01;
lastCommonTimePoint = 60;
i = 1;
AllParticles = alignedDatasetsStruct(i).AllParticles;
AllParticles(isnan(AllParticles))=0; %count off nuclei as having zero fluorescence
RepAbsTime = alignedDatasetsStruct(i).AbsTime;
Offset = size(AllParticles,2)-size(RepAbsTime,2)+1;
interpPoints = 500;
ReplicatesSpotsAccFluo = [];

for p = 1:size(AllParticles,1)
    particleFluo = AllParticles(p,Offset:end);
    interpTime = linspace(0,ceil(max(RepAbsTime)),interpPoints);
    interpFluo = interp1(RepAbsTime,particleFluo,interpTime);
    particleAccumulatedFluo = cumtrapz(interpTime,interpFluo);
    particleAccumulatedFluoDeg = [];
    particleAccumulatedFluoDeg(1) = 0;
    
    for t = 2:length(interpTime)
        mRNApreviousStep = particleAccumulatedFluoDeg(t-1);
        Production = particleAccumulatedFluo(t) - particleAccumulatedFluo(t-1);
        Degradation = -gamma * mRNApreviousStep;
        particleAccumulatedFluoDeg = [particleAccumulatedFluoDeg nansum([mRNApreviousStep;Degradation;Production])];       
    end
       
    finalSpotmRNA = particleAccumulatedFluoDeg(end);
    ReplicatesSpotsAccFluo(p) = finalSpotmRNA;
end
    
histogram(log10(ReplicatesSpotsAccFluo))
title(num2str(mean(ReplicatesSpotsAccFluo)./std(ReplicatesSpotsAccFluo)))
    