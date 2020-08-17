function OutputArray = CreateOffHomologDataArray(CompiledParticles,FrameInfo,Feature,lastCommonTime)

% Feature can be 'IntegratedFluo', 'MeanFluo' or 'Duration'
baselineAUs = 10; %the fluorescence assigned to undetected spots
minFranes = 3; % the shortest duration a spot can have to be considered
figure
%% * MEAN FLUORESCENCE *

if strcmpi(Feature,'MeanFluo')
    hold on
    usedParticles = [0];
    AllFirstParticlesMean = [];%nan(1,length(CompiledParticles));
    AllSecondParticlesMean = [];%nan(1,length(CompiledParticles));    
    for p1 = 1:length(CompiledParticles)
        p1
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        MeanFluo1 = mean(Fluo1);
        
        if length(Fluo1) > minFranes

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

            elseif HomParticle==0 || any(usedParticles==HomParticle) ||...
                CompiledParticles(HomParticle).HomologParticle ~= p1%~any(usedParticles==p1)% 
                FluoPair = [MeanFluo1 baselineAUs];
                FluoPair = FluoPair(randperm(length(FluoPair)));
                AllFirstParticlesMean = [AllFirstParticlesMean FluoPair(1)];
                AllSecondParticlesMean = [AllSecondParticlesMean FluoPair(2)];
                usedParticles = [usedParticles p1];         
            end

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
    LastCommonFrame = find(AbsTime<lastCommonTime,1,'last');

    
    hold on
    usedParticles = [0]; %to keep track of things
    AllFirstParticlesInt = [];%nan(1,length(CompiledParticles)); %to store data
    AllSecondParticlesInt = [];%nan(1,length(CompiledParticles));  
    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        Fluo1 = CompiledParticles(p1).Fluo;
        Frames1 = CompiledParticles(p1).Frame;
        lastFrame1 = find(Frames1 <= LastCommonFrame,1,'last');
        Fluo1 = Fluo1(1:lastFrame1);
        
        if length(Fluo1) > minFranes

            Frames1 = CompiledParticles(p1).Frame;
            Frames1 = Frames1(1:lastFrame1);
            %put background fluo in the frames where the particle doesn't exist
            ArrayToFill1 = zeros(1,length(Frames1(1):Frames1(end)))+baselineAUs;
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
                lastFrame2 = find(Frames2 <= LastCommonFrame,1,'last');
                Fluo2 = Fluo2(1:lastFrame2);
                Frames2 = Frames2(1:lastFrame2);
                
                %put background in the frames where the particle doesn't exist
                ArrayToFill2 = zeros(1,length(Frames2(1):Frames2(end)))+baselineAUs;
                %populate with fluorescence values
                FluoToFill2 = ArrayToFill2;
                TimeToFill2 = ArrayToFill2;
                FluoToFill2(Frames2) = Fluo2;
                TimeToFill2(Frames2) = AbsTime(Frames2);
                % now use trapezoidal integration on this trace
                IntFluo2 = trapz([TimeToFill2(1)-frameRate TimeToFill2 TimeToFill2(end)+frameRate]...
                    ,[0 FluoToFill2 0]);        
                % randomize who's allele #1 or #2 for plotting purposes
                FluoPair = [IntFluo1 IntFluo2];
                FluoPair = FluoPair(randperm(length(FluoPair)));
                AllFirstParticlesInt = [AllFirstParticlesInt FluoPair(1)];
                AllSecondParticlesInt = [AllSecondParticlesInt FluoPair(2)];
                usedParticles = [usedParticles p1 HomParticle];
            
            % now in case the particle doesn't have a homolog, assign a
            % 'fake' particle to it with a fluorescence equal to the
            % detection threshold
            elseif HomParticle == 0 || any(usedParticles==HomParticle)||...
                CompiledParticles(HomParticle).HomologParticle ~= p1 %if it doesn't have an allele in the same nucleus
                IntFluo2 = trapz([TimeToFill1(1)-frameRate TimeToFill1 TimeToFill1(end)+frameRate]...
                ,[0 ones(1,length(FluoToFill1))*baselineAUs 0]);
                FluoPair = [IntFluo1 IntFluo2];
                FluoPair = FluoPair(randperm(length(FluoPair)));
                AllFirstParticlesInt = [AllFirstParticlesInt FluoPair(1)];
                AllSecondParticlesInt = [AllSecondParticlesInt FluoPair(2)];
                usedParticles = [usedParticles p1];
            end
        end
    end
    
    OutputArray = [AllFirstParticlesInt;AllSecondParticlesInt];
    plot(AllFirstParticlesInt,AllSecondParticlesInt,'ro')
    hold on
    plot([0 max(OutputArray)*1.2],[0 max(OutputArray)*1.2],'k-')
    title('Integrated Particle Fluorescence Over Time')
end

%%  TIME ON
AbsTime = [FrameInfo.Time]./60;
lastFrame = AbsTime(end);

if strcmpi(Feature,'tOn')
    hold on
    usedParticles = [0];
    AllFirstParticlestON = [];
    AllSecondParticlestON = []; 
    
    for p1 = 1:length(CompiledParticles)
        HomParticle = CompiledParticles(p1).HomologParticle;
        tON1 = CompiledParticles(p1).FirstFrame;
        tON1 = AbsTime(tON1);
        
        %if they are mutual homologs and haven't been used already
        if HomParticle && CompiledParticles(HomParticle).HomologParticle == p1...
                && ~any(usedParticles==p1) && ~any(usedParticles==HomParticle) 
            tON2 = CompiledParticles(HomParticle).FirstFrame;
            tON2 = AbsTime(tON2);
            % randomize who's allele 1 and who's allele 2 for visualization purposes
            ValuePair = [tON1 tON2];
            ValuePair = ValuePair(randperm(length(ValuePair)));
            AllFirstParticlestON = [AllFirstParticlestON ValuePair(1)];
            AllSecondParticlestON = [AllSecondParticlestON ValuePair(2)];
            usedParticles = [usedParticles p1 HomParticle];

        elseif HomParticle==0 || any(usedParticles==HomParticle) ||...
            CompiledParticles(HomParticle).HomologParticle ~= p1%~any(usedParticles==p1)% 
            ValuePair = [tON1 lastFrame];
            ValuePair = ValuePair(randperm(length(ValuePair)));
            AllFirstParticlestON = [AllFirstParticlestON ValuePair(1)];
            AllSecondParticlestON = [AllSecondParticlestON ValuePair(2)];
            usedParticles = [usedParticles p1];         
        end

    end
        
    OutputArray = [AllFirstParticlestON;AllSecondParticlestON];
    plot(AllFirstParticlestON,AllSecondParticlestON,'ko')
    hold on
    plot([0 max(OutputArray)*1.2],[0 max(OutputArray)*1.2],'k-')
    title('time On')
        
end

    
end


