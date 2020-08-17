function OutputArray = CreateHomologDataArray(CompiledParticles,FrameInfo,Feature)

% Feature can be 'IntegratedFluo', 'MeanFluo' or 'Duration'

figure
%% * MEAN FLUORESCENCE *

if strcmpi(Feature,'MeanFluo')
    hold on
    usedParticles = [0];
    AllFirstParticlesMean = [];
    AllSecondParticlesMean = [];    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        MeanFluo1 = mean(Fluo1);
        
        %if they are mutual homologs and haven't been used already
        if HomParticle && CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) 
            Fluo2 = CompiledParticles(HomParticle).Fluo;
            MeanFluo2 = mean(Fluo2);
            % randomize who's allele 1 and who's allele 2 for visualization purposes
            FluoPair = [MeanFluo1 MeanFluo2];
            FluoPair = FluoPair(randperm(length(FluoPair)));
            AllFirstParticlesMean = [AllFirstParticlesMean FluoPair(1)];
            AllSecondParticlesMean = [AllSecondParticlesMean FluoPair(2)];
            usedParticles = [usedParticles p1 HomParticle];
        end
    end
    
    OutputArray = [AllFirstParticlesMean;AllSecondParticlesMean];
    plot(AllFirstParticlesMean,AllSecondParticlesMean,'ko')
    hold on
    plot([0 max(OutputArray)*1.2],[0 max(OutputArray)*1.2],'k-')
    title('Mean Particle Fluorescence Over Time')
end


%% * INTEGRATED FLUORESCENCE *
AbsTime = [FrameInfo.Time]./60;
frameRate = 1; %minutes per frame

if strcmpi(Feature,'IntegratedFluo')
    hold on
    usedParticles = [0]; %to keep track of thinga
    AllFirstParticlesInt = []; %to store data
    AllSecondParticlesInt = [];    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        Frames1 = CompiledParticles(p1).Frame;
        %put zeros in the frames where the particle doesn't exist
        ArrayToFill1 = zeros(1,length(Frames1(1):Frames1(end)));
        %populate with fluorescence and absolute time values
        FluoToFill1 = ArrayToFill1;
        TimeToFill1 = ArrayToFill1;
        FluoToFill1(Frames1) = Fluo1;
        TimeToFill1(Frames1) = AbsTime(Frames1);
        % now use trapezoidal integration on this trace
        IntFluo1 = trapz([TimeToFill1(1)-frameRate TimeToFill1 TimeToFill1(end)+frameRate]...
            ,[0 FluoToFill1 0]);
               
        %if they are mutual homologs and haven't been used already
        if HomParticle && CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) 
            Fluo2 = CompiledParticles(HomParticle).Fluo;
            Frames2 = CompiledParticles(HomParticle).Frame;
            %put zeros in the frames where the particle doesn't exist
            ArrayToFill2 = zeros(1,length(Frames2(1):Frames2(end)));
            %populate with fluorescence values
            FluoToFill2 = ArrayToFill2;
            TimeToFill2 = ArrayToFill2;
            FluoToFill2(Frames2) = Fluo2;
            TimeToFill2(Frames2) = AbsTime(Frames2);
            % now use trapezoidal integration on this trace
            IntFluo2 = trapz([TimeToFill2(1)-frameRate TimeToFill2 TimeToFill2(end)+frameRate]...
                ,[0 FluoToFill2 0]);        
            % randomize who's allele1 or 2
            FluoPair = [IntFluo1 IntFluo2];
            FluoPair = FluoPair(randperm(length(FluoPair)));
            AllFirstParticlesInt = [AllFirstParticlesInt FluoPair(1)];
            AllSecondParticlesInt = [AllSecondParticlesInt FluoPair(2)];
            usedParticles = [usedParticles p1 HomParticle];
        end
    end
    
    OutputArray = [AllFirstParticlesInt;AllSecondParticlesInt];
    plot(AllFirstParticlesInt,AllSecondParticlesInt,'ro')
    hold on
    plot([0 max(OutputArray)*1.2],[0 max(OutputArray)*1.2],'k-')
    title('Integrated Particle Fluorescence Over Time')
end

%% time to turn on

if strcmpi(Feature,'tON')
    hold on
    usedParticles = [0];
    AllFirstParticles_tON = [];
    AllSecondParticles_tON = [];    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        tON1 = CompiledParticles(p1).FirstFrame;
        
        %if they are mutual homologs and haven't been used already
        if HomParticle && CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) 
            tON2 = CompiledParticles(HomParticle).FirstFrame;
            % randomize who's allele 1 and who's allele 2 for plotting purposes
            tONPair = [tON1 tON2];
            tONPair = tONPair(randperm(length(tONPair)));
            AllFirstParticles_tON = [AllFirstParticles_tON tONPair(1)];
            AllSecondParticles_tON = [AllSecondParticles_tON tONPair(2)];
            usedParticles = [usedParticles p1 HomParticle];
        end
    end
    
    OutputArray = [AllFirstParticles_tON;AllSecondParticles_tON];
    plot(AllFirstParticles_tON,AllSecondParticles_tON,'ko')
    hold on
    plot([0 max(OutputArray)*1.2],[0 max(OutputArray)*1.2],'k-')
    title('tON')
end
















end
