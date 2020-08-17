function Particles = disapproveOverlappingParticles(Particles,Spots)

% Gets rid of particles occupying the exact same x,y position in the
% Particles.mat structure


CopyOfParticles = Particles;
ParticleMatrix = nan(length(Particles),length(Particles));

% make all posible pairwise comparisons between particles.
% store the pair distance in a square matrix where each particle occupies a
% given row and column.
for p1 = 1:length(Particles)
    particle1XPos = Particles(p1).xPos;
    particle1YPos = Particles(p1).yPos;
    for p2 = 1:length(Particles)
        particle2XPos = Particles(p2).xPos;
        particle2YPos = Particles(p2).yPos;        
        ParticleMatrix(p1,p2) = pdist([particle1XPos,particle1YPos ; particle2XPos,particle2YPos]);
    end
end

%imshow(ParticleMatrix,[])
% triu extracts the upper triangular part of a matrix
BadParticlesMatrix = triu(ParticleMatrix==0); %if a value outside the diagonal is 0 they are overlapping
DiagonalEntries = diag(diag(BadParticlesMatrix)); % we don't want to take the diagonal into account
BadParticlesMatrix(DiagonalEntries) = 0;
%imshow(BadParticlesMatrix)

% we want to find positions in the matrix that are 1
BadParticles = find(sum(BadParticlesMatrix,1)>0); %this returns the index of particles that had a distance of 0 to another one
% flag bad particles
for bp = BadParticles
    Particles(bp).Approved = -1;
end

% force approval of the rest
for p = 1:length(Particles)
    if Particles(p).Approved > -1
        Particles(p).Approved = 1;
    end
end

%fill z position if missing
for p = 1:length(Particles)
    %if isempty(Particles(p).zPos)
        Particles(p).zPos = ones(size(Particles(p).xPos));
   % end
end




end





