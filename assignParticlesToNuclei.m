function Nuclei = assignParticlesToNuclei (Nuclei,Mask,CompiledParticles)


% create a struct called 'Nuclei' to store the spots associated to each
% nucleus in the mask
Nuclei = regionprops(logical(Mask));
NucleusLabels = bwlabel(Mask);

% We'll assign particles to nuclei that exist in their position so, just ti be safe
% let's dilate the nuclei areas a little
SE = strel('disk',3);
DilatedNucleusLabels = imdilate(NucleusLabels,SE);


% sometimes it's easier to deal with row arrays instead of structs
NucleiXPos = zeros(1,length(Nuclei));
NucleiYPos = zeros(1,length(Nuclei));
for n = 1:length(Nuclei)
    NucleiYPos(n) = Nuclei(n).Centroid(2);
    NucleiXPos(n) = Nuclei(n).Centroid(1);
    Nuclei(n).AssociatedParticles1 = [];
    Nuclei(n).AssociatedParticles2 = [];
end
    
failcount = 0; %this keeps track of how many particles were left orphan
for p = 1:length(CompiledParticles)
    particleXPos = CompiledParticles(p).xPos;
    particleYPos = CompiledParticles(p).yPos;
    DistancesToEachNucleus = nan(1,length(Nuclei));
    
    % METHOD 1: assign particles based on proximity to nucleus centroid
    for n = 1:length(Nuclei)
        NucleusXPos = NucleiXPos(n);
        NucleusYPos = NucleiYPos(n);
        DistancesToEachNucleus(n) = sqrt((NucleusXPos-particleXPos)^2+(NucleusYPos-particleYPos)^2);
    end   
    [~,Idx] = min(DistancesToEachNucleus);
    Nuclei(Idx).AssociatedParticles1(length(Nuclei(Idx).AssociatedParticles1)+1) = p;
    
    %METHOD 2: assign particles based on overlap with nucleus area region
    if DilatedNucleusLabels(particleYPos,particleXPos) %if the particle overlaps with a nucleus area
        Nuclei(DilatedNucleusLabels(particleYPos,particleXPos)).AssociatedParticles2(length(Nuclei(Idx).AssociatedParticles2)+1) = p;
    else %find the closest nucleus area if the particle doesn't overlap with a nucleus
        %failcount = failcount+1;
        % select a square around the unassigned particle
        NearPixels = DilatedNucleusLabels(max(particleYPos-40,0):min(particleYPos+40,size(NucleusLabels,1)),...
            max(particleXPos-40,0):min(particleXPos+40,size(NucleusLabels,2)));
        [~,IDX] = bwdist(NearPixels); % for each pixel, calculate the distance to the closest nonzero pixel and its position
        NucleusID = NearPixels(IDX(40,40)); % 40,40 is the position of the particle in this square
        Nuclei(NucleusID).AssociatedParticles2(length(Nuclei(Idx).AssociatedParticles2)+1) = p;
    end       
    
end

% I want to compare both methods to make sure the results are identycal

% figure
% hold on
% scatter([CompiledParticles.xPos],[CompiledParticles.yPos],'r.')
% scatter(NucleiXPos,NucleiYPos,'bo')
% hold off
% 

%number of particles per nucleus
for n = 1:length(Nuclei)
    NumberOfParticles = length([Nuclei(n).AssociatedParticles1]);
    Nuclei(n).NParticles = NumberOfParticles;
end
% 
% save('Nuclei','Nuclei')
% save('DilatedNucleusLabels','DilatedNucleusLabels')

try
figure
scatter([Nuclei.AssociatedParticles1],[Nuclei.AssociatedParticles2],'ko')
catch
end

end