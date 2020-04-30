function [ErrorData,DimmestFluoData] =  FluoErrorDimmestSpots(Path,Struct)

% PATH is the path to the DynamicsResults folder
% STRUCT is a structure containing the compiled data from all replicates.
% It is generated by the CombineReplicatesFigures function

%for a given movie, what portion of the dimmest particles we take? it's going to be 1/N 
Ndimmestparticles = 4; % so =4 is the dimmest 0.25.

%for a dim particle, how many frames are we taking?
NdimmestFrames = 3;

ErrorData = [];
DimmestFluoData = [];


for p1 = 1:length(Struct) % first find what are the dimmest particles
    
    % load the CompiledParticles struct corresponding to this prefix
    Prefix = Struct(p1).Prefix;
    DynamicsResultsPath = [Path '/' Prefix];
    load([DynamicsResultsPath '/CompiledParticles.mat'],'CompiledParticles');
    if iscell(CompiledParticles) % in old versions of the lab code this could be a cell
        CompiledParticles = CompiledParticles{1};
    end
    
    %now find the particles with the lowest mean fluorescence
    particlesMeanFluo = nan(1,length(CompiledParticles));
    particleIndices = [1:length(CompiledParticles)];
    for p1 = 1:length(particlesMeanFluo)
       particlesMeanFluo(p1) = nanmean(CompiledParticles(p1).Fluo);
    end
    % this is to sort particles by fluorescence and extend that sorting to
    % their indices.
    A = [particlesMeanFluo;particleIndices]';
    B = sortrows(A,'ascend');
    % the indices of the dimmest quartile based on their mean fluo over time
    DimmestParticleIndices = B(ceil(1:length(CompiledParticles)/Ndimmestparticles),2); 

    % now use the indices to do the thing
    
    for p2 = 1:DimmestParticleIndices
        
        particleFluo = CompiledParticles(p2).Fluo;
        particleFluo = particleFluo(particleFluo>0); % sometimes there used to be spurious negative spot fluo fits
        dimmestFramesFluo = sort(particleFluo);
        dimmestFramesFluo = dimmestFramesFluo(1:min(NdimmestFrames,length(dimmestFramesFluo)));
        dimmestFramesFluo = mean(dimmestFramesFluo);
        particleError = CompiledParticles(p2).FluoError;
        
        ErrorData = [ErrorData particleError];
        DimmestFluoData = [DimmestFluoData dimmestFramesFluo];

    end
    
    
end

    
    
    