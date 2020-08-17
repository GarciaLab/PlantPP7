function PlotParticleScatter(CompiledParticles,Prefix,lastFrame,FrameInfo,Feature)

% makes a scatterplot of spots in the same nucleus and calculates total,
% correlated and intrinsic noise.
% FEATURE specifies which aspect/metric of gene expression we use for the
% noise calculations. it can be:
%'MeanFluo': the mean fluorescence over time. Inactive frames are not
%counted.
% 'IntegratedFluo' : area under the curve of the fluorescence time trace.
% 'Duration': number of frames the spot is detected
% 'TimeOn': frame of first detection.
AbsoluteTimePerFrame = [FrameInfo.Time]/60; %in minutes

if strcmpi(Feature,'MeanFluo')
    %% Make a figure of the scatter without normalizing the distributions
    F1 = figure;

    hold on
    usedParticles = [0];
    AllFirstParticlesMean = [];
    AllSecondParticlesMean = [];
    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        ParticleFrames1 = CompiledParticles(p1).Frame;
        % find the last frame that is larger than 'lastFrame'
        ParticleFrames1 = ParticleFrames1(ParticleFrames1 <= lastFrame);
        if ParticleFrames1

            LastParticleFrame1 = ParticleFrames1(end);    
            [~,LastFrameIdx1] = find(ParticleFrames1==LastParticleFrame1);
            Fluo1 = Fluo1(1:LastFrameIdx1);
            Fluo1 = mean(Fluo1);

            if CompiledParticles(HomParticle).HomologParticle == p1...
                    && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs

                Fluo2 = CompiledParticles(HomParticle).Fluo;
                ParticleFrames2 = CompiledParticles(HomParticle).Frame;
                ParticleFrames2 = ParticleFrames2(ParticleFrames2 <= lastFrame);
                if ParticleFrames2
                    LastParticleFrame2 = ParticleFrames2(end); 
                    [~,LastFrameIdx2] = find(ParticleFrames2==LastParticleFrame2);
                    Fluo2 = Fluo2(1:LastFrameIdx2);
                    Fluo2 = mean(Fluo2);
                    FluoPair = [Fluo1 Fluo2];
                    FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes                                                 
                    plot(FluoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')
                    usedParticles = [usedParticles p1 HomParticle];
                    AllFirstParticlesMean(end+1) = Fluo1;
                    AllSecondParticlesMean(end+1) = Fluo2;
                end
            elseif CompiledParticles(HomParticle).HomologParticle ~= p1...
                    && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs

                plot(Fluo1,0,'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')

            end
        end
    end
    hold off
    title('Mean fluorescence')

    
% Now the figure related to the noise calculations
% first, find the mean of the distribution (to normalize against)

figure

subplot(1,2,1)
hold on
usedParticles = [0];
MeanAll = mean([AllFirstParticlesMean AllSecondParticlesMean]);
RawSpots1 = [];
RawSpots2 = [];
Spots1Fluos = []; %store the value of the first spot of each nucleus, normalized by the mean of all spots
Spots2Fluos = []; %store the value of the second spot of each nucleus
counter = 1;
for p1 = 1:length(CompiledParticles)
    p1
    HomParticle = CompiledParticles(p1).HomologParticle;    
    Fluo1 = CompiledParticles(p1).Fluo;
    ParticleFrames1 = CompiledParticles(p1).Frame;
    % find the last frame that is larger than 'lastFrame'
    ParticleFrames1 = ParticleFrames1(ParticleFrames1 <= lastFrame);
    if ParticleFrames1
        LastParticleFrame1 = ParticleFrames1(end);    
        [~,LastFrameIdx1] = find(ParticleFrames1==LastParticleFrame1);
        Fluo1 = Fluo1(1:LastFrameIdx1);
        NormFluo1 = mean(Fluo1)/MeanAll;

        if CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs

            Fluo2 = CompiledParticles(HomParticle).Fluo;
            ParticleFrames2 = CompiledParticles(HomParticle).Frame;
            ParticleFrames2 = ParticleFrames2(ParticleFrames2 <= lastFrame);
            if ParticleFrames2
                LastParticleFrame2 = ParticleFrames2(end); 
                [~,LastFrameIdx2] = find(ParticleFrames2==LastParticleFrame2);
                Fluo2 = Fluo2(1:LastFrameIdx2);
                NormFluo2 = mean(Fluo2)/MeanAll;
                FluoPair = [NormFluo1 NormFluo2];
                FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes                                 
                plot(FluoPair(1),FluoPair(2),'bo','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')
                usedParticles = [usedParticles p1 HomParticle];

                RawSpots1(counter) = mean(CompiledParticles(p1).Fluo);
                RawSpots2(counter) = mean(CompiledParticles(HomParticle).Fluo);
                Spots1Fluos(counter) = 1-NormFluo1; %distance to the normalized mean
                Spots2Fluos(counter) = 1-NormFluo2;
                counter = counter+1;
            end
        elseif CompiledParticles(HomParticle).HomologParticle ~= p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs
            FluoPair = [NormFluo1 0];
            FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes                
            plot(FLuoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','none')

        end
    end
end

maxX = max(1+abs(Spots1Fluos))*1.2;
maxY = max(1+abs(Spots2Fluos))*1.2;
plot([0 maxX],[0 maxY],'k-')
hold off
xlim([0 maxX])
ylim([0 maxY])

NoiseTotal = 1/2 * (mean(Spots1Fluos.^2) + mean(Spots2Fluos.^2))
NoiseIntrinsic = 1/2 * mean((Spots1Fluos - Spots2Fluos).^2)
NoiseCorrelated = mean(Spots1Fluos.*Spots2Fluos)

title({Prefix,['\eta_{tot}^2 =' num2str(NoiseTotal) ' \eta_{int}^2 =' num2str(NoiseIntrinsic) ...
    ' \eta_{corr}^2 =' num2str(NoiseCorrelated)]})
xlabel('mean activity allele 1 (normalized to mean of all alleles)')
ylabel('mean activity allele 2 (normalized to mean of all alleles)')
    
    
    
    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    
    
    
    
elseif strcmpi(Feature,'IntegratedFluo')
    F2 = figure;
    hold on
    usedParticles = [0];
    AllFirstParticlesMean = [];
    AllSecondParticlesMean = [];
    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        ParticleFrames1 = CompiledParticles(p1).Frame;
        ParticleTimes1 = AbsoluteTimePerFrame(ParticleFrames1);
        % find the last frame that is larger than 'lastFrame'
        ParticleFrames1 = ParticleFrames1(ParticleFrames1 <= lastFrame);
        if length(ParticleFrames1)

            LastParticleFrame1 = ParticleFrames1(end);    
            [~,LastFrameIdx1] = find(ParticleFrames1==LastParticleFrame1);
            ParticleTimes1 = ParticleTimes1(1:LastFrameIdx1);
            Fluo1 = Fluo1(1:LastFrameIdx1);
            Fluo1 = trapz([ParticleTimes1(1)-1 ParticleTimes1 ParticleTimes1(1)+1],...
                [ 0 Fluo1 0]);

            if CompiledParticles(HomParticle).HomologParticle == p1...
                    && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs

                Fluo2 = CompiledParticles(HomParticle).Fluo;
                ParticleFrames2 = CompiledParticles(HomParticle).Frame;
                ParticleFrames2 = ParticleFrames2(ParticleFrames2 <= lastFrame);
                ParticleTimes2 = AbsoluteTimePerFrame(ParticleFrames2);
                
                if ParticleFrames2
                    LastParticleFrame2 = ParticleFrames2(end); 
                    [~,LastFrameIdx2] = find(ParticleFrames2==LastParticleFrame2);
                    Fluo2 = Fluo2(1:LastFrameIdx2);
                    ParticleTimes2 = ParticleTimes2(1:LastFrameIdx2);
                    Fluo2 = trapz([ParticleTimes2(1)-1 ParticleTimes2 ParticleTimes2(end)+1],...
                        [0 Fluo2 0]);
                    FluoPair = [Fluo1 Fluo2];
                    FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes               
                    plot(FluoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')
                    usedParticles = [usedParticles p1 HomParticle];
                    AllFirstParticlesMean(end+1) = Fluo1;
                    AllSecondParticlesMean(end+1) = Fluo2;
                end
            elseif CompiledParticles(HomParticle).HomologParticle ~= p1...
                    && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs
                FluoPair = [Fluo1 0];
                FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes
                plot(FluoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')

            end
        end
    end
    hold off
    title('Integrated Fluorescence')

% now the figure related to the noise calculations
% first, find the mean of the distribution (to normalize against)

figure
%subplot(1,2,1)
hold on
usedParticles = [0];
MeanAll = mean([AllFirstParticlesMean AllSecondParticlesMean]);
RawSpots1 = [];
RawSpots2 = [];
Spots1Fluos = []; %store the value of the first spot of each nucleus, normalized by the mean of all spots
Spots2Fluos = []; %store the value of the second spot of each nucleus
counter = 1;
for p1 = 1:length(CompiledParticles)
    p1
    HomParticle = CompiledParticles(p1).HomologParticle;    
    Fluo1 = CompiledParticles(p1).Fluo;
    ParticleFrames1 = CompiledParticles(p1).Frame;
    % find the last frame that is larger than 'lastFrame'
    ParticleFrames1 = ParticleFrames1(ParticleFrames1 <= lastFrame);
    ParticleTimes1 = AbsoluteTimePerFrame(ParticleFrames1);
    if ParticleFrames1
        LastParticleFrame1 = ParticleFrames1(end);    
        [~,LastFrameIdx1] = find(ParticleFrames1==LastParticleFrame1);
        Fluo1 = Fluo1(1:LastFrameIdx1);
        Fluo1 = trapz([ParticleTimes1(1)-1 ParticleTimes1 ParticleTimes1(1)+1],...
                [ 0 Fluo1 0]);
        NormFluo1 = Fluo1/MeanAll;

        if CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs

            Fluo2 = CompiledParticles(HomParticle).Fluo;
            ParticleFrames2 = CompiledParticles(HomParticle).Frame;
            ParticleFrames2 = ParticleFrames2(ParticleFrames2 <= lastFrame);
            ParticleTimes2 = AbsoluteTimePerFrame(ParticleFrames2);
            if ParticleFrames2
                LastParticleFrame2 = ParticleFrames2(end); 
                [~,LastFrameIdx2] = find(ParticleFrames2==LastParticleFrame2);
                Fluo2 = Fluo2(1:LastFrameIdx2);
                ParticleTimes2 = ParticleTimes2(1:LastFrameIdx2);
                Fluo2 = trapz([ParticleTimes2(1)-1 ParticleTimes2 ParticleTimes2(end)+1],...
                        [0 Fluo2 0]);              
                NormFluo2 = Fluo2/MeanAll;
                FluoPair = [NormFluo1 NormFluo2];
                FluoPair = FluoPair(randperm(length(FluoPair))); %randomize for visualization purposes
                plot(FluoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','none')
                usedParticles = [usedParticles p1 HomParticle];
%                 RawSpots1(counter) = Fluo1;
%                 RawSpots2(counter) = Fluo2;
                Spots1Fluos(counter) = 1-NormFluo1; %distance to the normalized mean
                Spots2Fluos(counter) = 1-NormFluo2;
                counter = counter+1;
            end
        elseif CompiledParticles(HomParticle).HomologParticle ~= p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) %if they are mutual homologs
            FluoPair = [NormFluo1,0];
            FluoPair = FluoPair(randperm(length(FluoPair)));
            %plot(FluoPair(1),FluoPair(2),'o','MarkerSize',12,'MarkerFaceColor','r','MarkerEdgeColor','none')

        end
    end
end
maxX = max(1+abs(Spots1Fluos))*1.2;
maxY = max(1+abs(Spots2Fluos))*1.2;
plot([0 maxX],[0 maxY],'k-')
hold off
xlim([0 maxX])
ylim([0 maxY])

NoiseTotal = 1/2 * (mean(Spots1Fluos.^2) + mean(Spots2Fluos.^2))
NoiseIntrinsic = 1/2 * mean((Spots1Fluos - Spots2Fluos).^2)
NoiseCorrelated = mean(Spots1Fluos.*Spots2Fluos)

title({Prefix,['\eta_{tot}^2 =' num2str(NoiseTotal) ' \eta_{int}^2 =' num2str(NoiseIntrinsic) ...
    ' \eta_{corr}^2 =' num2str(NoiseCorrelated)]})
xlabel('integrated fluorescence allele 1 (normalized to mean of all alleles)')
ylabel('integrated fluorescence allele 2 (normalized to mean of all alleles)')





    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    
    
    
    
elseif strcmpi(Feature,'somefeature')






    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    
    
    
        
elseif strcmpi(Feature,'somefeature')
    
 




    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    % *\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*\/*
    
    
    
           
end
    
%% bootstrap the errors to get errobars
nSamples = 10000; %how many times we're bootstrapping
%these are the function we are bootstrapping, the formulas for each noise component
TotObsP = @(x,y) 1/2 * (mean(x.^2) + mean(y.^2)); 
IntObsP = @(x,y) 1/2 * mean((x-y).^2);
CorrObsP = @(x,y) mean(x.*y);

% bootstrapping the observed total noise
[bootTotObs,~] = bootstrp(nSamples,TotObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanTotalNoise = mean(bootTotObs);
bootstrappedStdTotalNoise = std(bootTotObs);

% bootstrapping the observed intrinsic noise
[bootIntObs,~] = bootstrp(nSamples,IntObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanIntrinsicNoise = mean(bootIntObs);
bootstrappedStdIntrinsicNoise = std(bootIntObs);

% bootstrapping the observed correlated noise
[bootCorrObs,~] = bootstrp(nSamples,CorrObsP,Spots1Fluos,Spots2Fluos); 
bootstrappedMeanCorrNoise = mean(bootCorrObs);
bootstrappedStdCorrNoise = std(bootCorrObs);

figure
%subplot(1,2,2)
hold on
errorbar(1,bootstrappedMeanTotalNoise,bootstrappedStdTotalNoise,'bo','LineWidth',3,...
    'CapSize',0)
errorbar(2,bootstrappedMeanIntrinsicNoise,bootstrappedStdIntrinsicNoise,'ro','LineWidth',3,...
    'CapSize',0)
errorbar(3,bootstrappedMeanCorrNoise,bootstrappedStdCorrNoise,'go','LineWidth',3,...
    'CapSize',0)
hold off
ylim([0 1])
xlim([0.8 3.2])
legend('\eta_{tot}^2','\eta_{int}^2','\eta_{corr}^2')




end
