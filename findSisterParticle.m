function CompiledParticles = findSisterParticle(CompiledParticles)

distanceThresh = 60; %in pixels, minimum distance for two particles to be considered as candidates to be in the same nucleus

for p1 = 1:length(CompiledParticles)
    
    FirstPartFrames = CompiledParticles(p1).Frame;
    ClosestParticlePerFrame = nan(1,length(FirstPartFrames));
    
    for f1 = 1:length(FirstPartFrames)
        actualFrame = FirstPartFrames(f1);
        Xpos1 = CompiledParticles(p1).xPos(f1);
        Ypos1 = CompiledParticles(p1).yPos(f1);
        DistanceToRestThisFrame = nan(1,length(CompiledParticles));
        
        
        for p2 = 1:length(CompiledParticles)
            if p2 ~= p1
            SecPartFrames = CompiledParticles(p2).Frame;
            [~,frIdx] = find(SecPartFrames == actualFrame);
            
            if frIdx
                Xpos2 = CompiledParticles(p2).xPos(frIdx);
                Ypos2 = CompiledParticles(p2).yPos(frIdx);

                DistanceP1P2thisFrame = sqrt((Xpos2-Xpos1)^2+(Ypos2-Ypos1)^2);
                
                if DistanceP1P2thisFrame < distanceThresh
                    DistanceToRestThisFrame(p2) = DistanceP1P2thisFrame;
                end
                
            end
            end
        end
        
        if ~isnan(min(DistanceToRestThisFrame))
            [~,ClosestParticlePerFrame(f1)] = min(DistanceToRestThisFrame);
        else
            ClosestParticlePerFrame(f1) = 0;
        end
    end
    
    CompiledParticles(p1).HomologParticle = mode(ClosestParticlePerFrame);

end



end
