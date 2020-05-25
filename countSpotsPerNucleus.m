function ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold)

% creates a 1xn array where n is the number of nuclei. It's populated with
% the number of spots in that nucleus

ParticlesPerNucleus = nan(1,length(Nuclei));

nucleicounter = 1;

for n = 17:length(Nuclei)
    NucleusParticles = Nuclei(n).AssociatedParticles1;
    NucleusArea = Nuclei(n).Area;
    
    if NucleusArea < AreaThreshold
        if ~isempty(NucleusParticles)
            ParticlesPerNucleus(nucleicounter) = length(NucleusParticles);
        else
            ParticlesPerNucleus(nucleicounter) = 0;
        end
        nucleicounter = nucleicounter+1;
    end
end
% figure to visualize output
%histogram(ParticlesPerNucleus)
end

