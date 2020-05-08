function ParticlesPerNucleus = countSpotsPerNucleus(Nuclei,AreaThreshold)

% creates a 1xn array where n is the number of nuclei. It's populated with
% the number of spots in that nucleus

ParticlesPerNucleus = nan(1,length(Nuclei));

nucleicounter = 0;

for n = 1:length(Nuclei)
    NucleusParticles = Nuclei(n).AssociatedParticles2;
    NucleusArea = Nuclei(n).Area;
    
    if NucleusArea < AreaThreshold
        nucleicounter = nucleicounter+1;
        if isempty(NucleusParticles)
            ParticlesPerNucleus(n) = length(NucleusParticles);
        else
            ParticlesPerNucleus(n) = 0;
        end
    end
end
