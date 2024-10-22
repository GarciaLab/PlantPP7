function particleVariability(Struct,Prefixes,DynamicsResultsPath,Feature,peakTime,PolPerAU,gamma)


% Feature can be 'derivative', 'fluorescence' or...

if strcmpi(Feature,'derivative')
    
    % for each particle, take the first derivative of its fluorescence over time
    % and divide each dFluo/dt(t) by the Fluo(t). Then get the frequency of these
    % values for all particles (so histograms of normalized
    % fluorescence derivative). We then take the mean and error of these
    % frequencies (or probabilities) to get a histogram with error.
    
    figure
    hold on
    BinWidth = .05;
    Edges = [-.5:BinWidth:.5]; %for fluo derivative

    HistogramCounts = nan(1000,length(Edges)-1);

    counter = 1;
    for rep = 1:length(Prefixes)
        prefixName = Prefixes{rep};
        prefixResultsFolder = [DynamicsResultsPath '/' prefixName];
        load([prefixResultsFolder '/' 'CompiledParticles.mat'],'CompiledParticles')
        if iscell(CompiledParticles) %in some versions of the analysis pipeline the struct is inside of a cell
            CompiledParticles = CompiledParticles{1}; %in that case we want the actual struct.
        end

        % now loop over this replicate's particles to calculate things
        minFrames = 10;% number of frames of the shortest particle analyzed
        for p = 1:length(CompiledParticles)
            particleFluo = smooth(CompiledParticles(p).Fluo,3);
            fluoDerivative = diff(particleFluo);
            % take the mean of each pair of consecutive timepoints:
            % we can normalize the derivative at time t by the mean between t-1
            % and t+1
            fluoConsecutiveMeans =  0.5 * (particleFluo(1:end-1) + particleFluo(2:end));

            if length(particleFluo)> minFrames
                [N,~] = histcounts(fluoDerivative(5:end)./fluoConsecutiveMeans(5:end),Edges,...
                    'Normalization','probability');
                HistogramCounts(counter,:) = N;
                counter = counter+1;
            end
        end
    end
    
    plotEdges = Edges(1:end-1);% + BinWidth;
    errorbar(plotEdges,nanmean(HistogramCounts,1),nanstd(HistogramCounts,1)./sqrt(counter),...
        'k','LineWidth',2,'CapSize',0)
    hold off
    xlabel('dFluo/dt (normalized by Fluo(t))')
    ylabel('frequency in time trace')
    xlim([plotEdges(1) plotEdges(end-1)])
    
elseif strcmpi(Feature,'fluorescence')
    % for each particle, take its fluorescence over time and divide by the
    % mean of its fluorescence over time. Then get the frequency of these
    % fluo(t) values for all particles (so histograms of normalized
    % fluorescence over time). We then take the mean and error of these
    % frequencies (or probabilities) to get a histogram with error.
   
    figure
    hold on
    BinWidth = .1;
    Edges = [-0.1:BinWidth:2.1]; 
    counter = 1;
    for rep = 1:length(Prefixes)
        prefixName = Prefixes{rep};
        prefixResultsFolder = [DynamicsResultsPath '/' prefixName];
        load([prefixResultsFolder '/' 'CompiledParticles.mat'],'CompiledParticles')
        if iscell(CompiledParticles) %in some versions of the analysis pipeline the struct is inside of a cell
            CompiledParticles = CompiledParticles{1}; %in that case we want the actual struct.
        end

        % now loop over this replicate's particles to calculate things
        minFrames = 10;% number of frames of the shortest particle analyzed
        for p = 1:length(CompiledParticles)
            particleFluo = smooth(CompiledParticles(p).Fluo,3);
            particleMeanFluo = mean(particleFluo);

            if length(particleFluo)> minFrames
                [N,~] = histcounts(particleFluo(5:end)./particleMeanFluo,Edges,...
                    'Normalization','probability');
                HistogramCounts(counter,:) = N;
                ForCVcalculation(counter) = std(particleFluo(5:end))./particleMeanFluo;
                counter = counter+1;
            end
            
        end
    end
    CV = mean(ForCVcalculation);
    plotEdges = Edges(1:end-1);% + BinWidth;
    errorbar(plotEdges,nanmean(HistogramCounts,1),nanstd(HistogramCounts,1)./sqrt(counter),...
        'k','LineWidth',2,'CapSize',0)
    plot([1-CV 1+CV],[0.25 0.25],'k')
    plot(1,0.25,'ko')
    hold off
    xlabel('fluo(t) (normalized by mean trace fluo)')
    ylabel('frequency in time trace')
    xlim([plotEdges(1) plotEdges(end-1)])
    
    
elseif strcmpi(Feature,'AllParticles')
    
    AllParticlesPeakArray = nan(sum(cell2mat(cellfun(@(x) size(x, 1),{Struct.AllParticles},'UniformOutput',false))),...
        4); 
    % we'll store things in this array of px3 where p is the total number in all 
    % the replicates which in turn is equivalent to sum of alignedDatasetsStruct.AllParticles 
    % matrices along the first dimension.
    
    firstrow = 1;
    for rep = 1:length(Prefixes)
        
        % first find which column correspond to the timepoint we want to
        % plot (specified by peakTime)
        timeDiff = Struct(rep).AbsTime - peakTime;
        ClosestFrame = find(abs(timeDiff) == min(abs(timeDiff)));
        
        AllRepParticles = Struct(rep).AllParticles;
        AllParticlesPeakArray(firstrow:firstrow+size(AllRepParticles,1)-1,:) = (AllRepParticles(:,ClosestFrame-2:ClosestFrame+1));
        firstrow = firstrow + size(AllRepParticles,1);
    end
    % convert 0s to nans to show only actively transcribing loci
    AllParticlesPeakArray(AllParticlesPeakArray<=0) = nan;
    meanFluoAcrossThese4Frames = nanmean(AllParticlesPeakArray(:));
    NormAllParticles = AllParticlesPeakArray./meanFluoAcrossThese4Frames;
    
    %Edges = [0:0.05:8];
    Edges = logspace(-1.6,1.1,20);
    %histogram(NormAllParticles,Edges,'Normalization','probability')
    [N,~] = histcounts(NormAllParticles,Edges,'Normalization','probability');
    hold on
    histogram(NormAllParticles,Edges,'Normalization','probability','EdgeColor','none')
    CV = mean(nanstd(NormAllParticles));
    plot([1-CV/2 1+CV/2],[0.1 0.1],'k-')
    plot(1,0.1,'ko')
    %plot(Edges(1:end-1)+(Edges(2:end)/10),N,'b-','LineWidth',2);
    set(gca,'xscale','log')
    hold off
   

elseif strcmpi(Feature,'AllParticlesMean')
    % calculate the distribution of mean particle fluorescence over time
    % we'll put all the fluo traces from all replicates in this struct
    % frames with zero fluorescence won't be counted (will be made nans)
    AllParticles =[];
    for rep = 1:length(Struct)
        AllReplicateParticles = Struct(rep).AllParticles;
        
        % we will count only the frames from first to last detection.
        % missing frames will be counted as nans
        for p = 1:size(AllReplicateParticles,1)
            particleTrace = AllReplicateParticles(p,:);
            firstNonZeroFrame = find(particleTrace,1);
            lastNonZeroFrame = find(particleTrace,1,'last');
            % make the frames prior to detection nans
            AllReplicateParticles(p,1:firstNonZeroFrame-1) = nan;
            % make the frames after detection nans
            AllReplicateParticles(p,lastNonZeroFrame+1:end) = nan;
            % make the frames in between that are undetected nans
            AllReplicateParticles(AllReplicateParticles==0) = nan;
        end
        
        MeanOfEachParticle = nanmean(AllReplicateParticles,2); % mean across frames (columns)
        MeanOfEachParticle = MeanOfEachParticle(~isnan(MeanOfEachParticle));
        AllParticles = [AllParticles;MeanOfEachParticle];
    end
    
    figure
    % we count only the frames with signal
    Mean = nanmean(AllParticles);
    NormAllParticles = AllParticles/Mean;
    SD = nanstd(NormAllParticles);
    CV = SD;
    histogram(NormAllParticles,logspace(-1.3,1.8,19),'Normalization','Probability',...
        'EdgeColor','k','DisplayStyle','Stairs','LineWidth',2);
    hold on
    plot([Mean-SD Mean+SD], [0.05 0.05],'k-')
    plot(Mean,0.05,'ko')
    hold off
    set(gca,'xscale','log')
    legend(['CV = ' num2str(CV)])
    
    figure
    AbsMean = nanmean(AllParticles.*PolPerAU);
    %AbsNormAllParticles = (AllParticles.*PolPerAU)/AbsMean;
    SD = nanstd(AllParticles.*PolPerAU);
    CV = SD/AbsMean;
    Hist = histogram(log10(AllParticles.*PolPerAU),18,'Normalization','Probability',...
        'EdgeColor','k','DisplayStyle','Stairs','LineWidth',2);
    hold on
    plot(log10([AbsMean-SD/2 AbsMean+SD/2]), [0.1 0.1],'k-')
    plot(log10(AbsMean),0.1,'ko')
    
    % gaussian fit to be fancier
%     % z(1), z(2), z(3) and z(4) (amplitude, center, spread and baseline)
%     DataX = mean([Hist.BinEdges(1:end-1);Hist.BinEdges(2:end)]);
%     DataY = Hist.Values;
%     Gaussfun = @(z)z(1) * exp(-(DataX-z(2)).^2./(2*(z(3)^2))) + z(4) - DataY; 
%     z0 = [0.02,log10(AbsMean),log10(AbsMean),0];
%     lb = [0.01,log10(AbsMean)*0.7,0.1,0]; % define lower bounds for each parameter
%     ub = [1,log10(AbsMean)*1.1,5,0]; % define upper bounds for each parameter, note that I'm forcing the baseline to be 0
%     Guess = lsqnonlin(Gaussfun,z0,lb,ub); %fit by minizing the subtraction between actual Y data and a gaussian
%     FitGauss = Guess(1)*exp(-(DataX-Guess(2)).^2./(2*(Guess(3)^2))) + Guess(4);
%     plot(DataX,FitGauss,'b','LineWidth',2)

    % spline fit
    splineRes = 0.1;
    xE = [Hist.BinEdges(1)-Hist.BinWidth/2 Hist.BinEdges(1:end-1)+Hist.BinWidth/2 Hist.BinEdges(end)];
    yE = [0 Hist.Values 0];
    xxE = xE(1):splineRes:xE(end);
    sE = spline(xE,yE,xxE);
    pE = pchip(xE,yE,xxE);
    plot(xxE,sE,'b','LineWidth',2)

    hold off
    %set(gca,'xscale','log')
    legend(['CV = ' num2str(CV)])
    title('mean RNAP occupancy across time')
    xlabel('log_{10}RNAP')
    ylabel('frequency')
    
    
    
figure
Palette = brewermap(length(Prefixes),'Dark2');
hold on
for rep = 1:length(Prefixes)
    Color = Palette(rep,:);
    ThisRepParticles = Struct(rep).AllParticles./Mean;
    %ParticlesMeansOverTime = mean(ThisRepParticles,2);
    %ThisRepParticles = ThisRepParticles./ParticlesMeansOverTime;
    ThisRepAbsTime = Struct(rep).AbsTime;
    
    for p = 1:size(ThisRepParticles,1)
        particleTrace = ThisRepParticles(p,:);
        if sum(~isnan(particleTrace))
            firstNonZeroFrame = find(particleTrace,1);
            Y = particleTrace(firstNonZeroFrame:end);
            X = ThisRepAbsTime(1:length(Y))
            X = X-X(1);
            plot(X,Y,'Color',Color)
        end
    end    
    %ThisRepParticles = ThisRepParticles(:,size(ThisRepParticles,2)-size(ThisRepAbsTime,2)+1:end)
    %ThisRepAbsTime = ThisRepAbsTime(size(ThisRepAbsTime,1)+1-size(ThisRepParticles,2):end);
    %plot(ThisRepAbsTime,ThisRepParticles,'Color',Color)
end
set(gca,'yscale','log')
%ylim([0.0025 20])
hold off


elseif strcmpi(Feature,'Integral')
    
    % not all movies last for the same ammount of time. For integration we have
    % to take that into account.
    LastFramePerPrefix = [];
    for p = 1:length(Struct)
        lastTimePoint = Struct(p).AbsTime;
        LastFramePerPrefix = [LastFramePerPrefix lastTimePoint(end)];
    end
    lastCommonTimePoint = min(LastFramePerPrefix);
    minAU = 0.1; %detection threshold AUs
    
    % for each replicate, find the frame that is closest to the last common
    % timepoint and get the integrated mRNA per spot up until that frame
    AllIntegratedFluo = [];
    for p = 1:length(Struct)
         integrationFrame = find(abs(Struct(p).AbsTime - lastCommonTimePoint) ...
             == min(abs(Struct(p).AbsTime - lastCommonTimePoint)));
         integratedFluo = Struct(p).IntegralSoFarWithOff;
         integratedFluo = integratedFluo(:,integrationFrame);
         undetected = integratedFluo(:,end)==0;
         integratedFluo(undetected) = lastCommonTimePoint*minAU;
         AllIntegratedFluo = [AllIntegratedFluo;integratedFluo];
    end
    
    
    %now do the histogram and all that
    figure
    Mean = nanmean(AllIntegratedFluo);
    %AbsNormAllParticles = (AllParticles.*PolPerAU)/AbsMean;
    SD = nanstd(AllIntegratedFluo);
    CV = SD/Mean;
    LogBinEdges = logspace(0.5,5,19);
    Hist = histogram(AllIntegratedFluo,LogBinEdges,'Normalization','Probability',...
        'EdgeColor','k','DisplayStyle','Stairs','LineWidth',2);
    hold on
    plot([Mean-(SD/2) Mean+(SD/2)], [0.1 0.1],'r-')
    plot(Mean,0.1,'bo')
%     % gaussian fit to be fancier
%     % z(1), z(2), z(3) and z(4) (amplitude, center, spread and baseline)
%     DataX = [Hist.BinEdges(1)-Hist.BinWidth/2 mean([Hist.BinEdges(1:end-1);Hist.BinEdges(2:end)]) Hist.BinEdges(end)+Hist.BinWidth/2];
%     DataY = [0 Hist.Values 0];
%     Gaussfun = @(z)z(1) * exp(-(DataX-z(2)).^2./(2*(z(3)^2))) + z(4) - DataY; 
%     z0 = [0.02,log10(Mean),log10(Mean),0];
%     lb = [0.01,log10(Mean)*0.2,0.1,0]; % define lower bounds for each parameter
%     ub = [1,log10(Mean)*1.5,5,0]; % define upper bounds for each parameter, note that I'm forcing the baseline to be 0
%     Guess = lsqnonlin(Gaussfun,z0,lb,ub); %fit by minizing the subtraction between actual Y data and a gaussian
%     FitGauss = Guess(1)*exp(-(DataX-Guess(2)).^2./(2*(Guess(3)^2))) + Guess(4);
%     plot(DataX,FitGauss,'b','LineWidth',2)
    hold off
    %set(gca,'xscale','log')
    legend(['CV = ' num2str(CV)])
    title('integrated spot fluorescence')
    xlabel('log_{10} integrated fluorescence')
    ylabel('frequency')
    set(gca,'XScale','log')



    
    
elseif strcmpi(Feature,'IntegralWithDegradation')
    
    LastFramePerPrefix = [];
    for pfx = 1:length(Struct)
        lastTimePoint = Struct(pfx).AbsTime;
        LastFramePerPrefix = [LastFramePerPrefix lastTimePoint(end)];
    end
    lastCommonTimePoint = min(LastFramePerPrefix);
    AllReplicatesSpotsAccFluo = [];
    
    
    for pfx = 1:length(Struct)
        
        AllParticles = Struct(pfx).AllParticles;
        AllParticles(isnan(AllParticles))=0; %count off nuclei as having zero fluorescence
        RepAbsTime = Struct(pfx).AbsTime;
        Offset = size(AllParticles,2)-size(RepAbsTime,2)+1;
        interpPoints = 500;
        ThisReplicateSpotsAccFluo = [];

        for p = 1:size(AllParticles,1)
            particleFluo = AllParticles(p,Offset:end);
            interpTime = linspace(0,ceil(max(RepAbsTime)),interpPoints);
            interpFluo = interp1(RepAbsTime,particleFluo,interpTime);
            particleAccumulatedFluo = cumtrapz(interpTime,interpFluo);
            particleAccumulatedFluoDeg = [];
            particleAccumulatedFluoDeg(1) = 0;

            for t = 2:length(interpTime)
                mRNApreviousStep = particleAccumulatedFluoDeg(t-1);
                Production = particleAccumulatedFluo(t) - particleAccumulatedFluo(t-1);
                Degradation = -gamma * mRNApreviousStep;
                particleAccumulatedFluoDeg = [particleAccumulatedFluoDeg nansum([mRNApreviousStep;Degradation;Production])];       
            end

            finalSpotmRNA = particleAccumulatedFluoDeg(end);
            ThisReplicateSpotsAccFluo(p) = finalSpotmRNA;
        end
        AllReplicatesSpotsAccFluo = [AllReplicatesSpotsAccFluo ThisReplicateSpotsAccFluo];
    end
    
    MinPerTimeStep = (interpTime(3)-interpTime(2)); %in min
    AbsGamma = gamma * (1/MinPerTimeStep); % in min-1
    halfLife = log(2)/AbsGamma; %in minutes
    
    
    
    histogram(log10(AllReplicatesSpotsAccFluo),20,'Normalization','Probability');
    
    title(['CV= ' num2str(std(AllReplicatesSpotsAccFluo)./mean(AllReplicatesSpotsAccFluo)) ...
        '  %Off= ' num2str(sum((AllReplicatesSpotsAccFluo==0))./length(AllReplicatesSpotsAccFluo)) ...
        '  Half-life(min)= ' num2str(halfLife)])
end


    
    

end

 
    
   
    
    



    
    


    







