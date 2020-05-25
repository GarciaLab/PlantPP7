function [NucleiPerFrame, CompiledParticles] = ...
    assignParticlesToEllipses (Ellipses,CompiledParticles,FrameInfo)


for frame = 1:length(FrameInfo)

    % create the Nuclei struct for output
    Nuclei = struct('AssociatedParticles1',[],'NParticles',[],'Area',[]);
    ThisFrameParticles = struct('xPos',[],'yPos',[],'OriginalParticle',[]);
    
    ThisFrameEllipses = Ellipses{frame};
    
    for e = 1:length(ThisFrameEllipses)
        Nuclei(e).AssociatedParticles1 = [];
        Nuclei(e).Area = 0;
    end
    
    NucleiXPos = ThisFrameEllipses(:,1);
    NucleiYPos = ThisFrameEllipses(:,2);
    
    particleCounter = 1;
    % grab the particles that exist in this frame
    for cp = 1:length(CompiledParticles)
        particleFrames = CompiledParticles(cp).Frame;
        [~,idx] = find(particleFrames == frame); %index of the frame in the particle's life
        if idx
            ThisFrameParticles(particleCounter).xPos = CompiledParticles(cp).xPos(idx);
            ThisFrameParticles(particleCounter).yPos = CompiledParticles(cp).yPos(idx);
            ThisFrameParticles(particleCounter).OriginalParticle = cp;
            particleCounter = particleCounter+1;
        end
    end


    failcount = 0; %this keeps track of how many particles were left orphan
    if ~isempty(ThisFrameParticles(1).xPos)
    for p = 1:length(ThisFrameParticles)
        particleXPos = ThisFrameParticles(p).xPos;
        particleYPos = ThisFrameParticles(p).yPos;
        DistancesToEachNucleus = nan(1,length(NucleiXPos));

        % METHOD 1: assign particles based on proximity to nucleus centroid
        for n = 1:length(NucleiXPos)
            NucleusXPos = NucleiXPos(n);
            NucleusYPos = NucleiYPos(n);
            DistancesToEachNucleus(n) = sqrt((NucleusXPos-particleXPos)^2+(NucleusYPos-particleYPos)^2);
        end   
        [~,Idx] = min(DistancesToEachNucleus);
        Nuclei(Idx).AssociatedParticles1 = [Nuclei(Idx).AssociatedParticles1 p];

        %Add the ellipse info to the CompiledParticles struct, we might want it
        %later on
        OriginalParticleIdx = ThisFrameParticles(p).OriginalParticle;
        CompiledParticles(OriginalParticleIdx).AssociatedEllipse(frame) = Idx;

    end
    end

    %number of particles per nucleus
    for n = 1:length(Nuclei)
        NumberOfParticles = length([Nuclei(n).AssociatedParticles1]);
        Nuclei(n).NParticles = NumberOfParticles;
    end
    
    NucleiPerFrame(frame).Nuclei = Nuclei;
end





end