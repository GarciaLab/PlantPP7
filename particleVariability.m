function particleVariability(Struct,Prefixes,DynamicsResultsPath,Feature,peakTime)


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
                counter = counter+1;
            end
            
        end
    end
    
    plotEdges = Edges(1:end-1);% + BinWidth;
    errorbar(plotEdges,nanmean(HistogramCounts,1),nanstd(HistogramCounts,1)./sqrt(counter),...
        'k','LineWidth',2,'CapSize',0)
    hold off
    xlabel('fluo(t) (normalized by mean trace fluo)')
    ylabel('frequency in time trace')
    xlim([plotEdges(1) plotEdges(end-1)])
    

    
    
    
elseif strcmpi(Feature,'AllParticles')
    
    AllParticlesArray = nan(sum(cell2mat(cellfun(@(x) size(x, 1),{Struct.AllParticles},'UniformOutput',false))),...
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
        AllParticlesArray(firstrow:firstrow+size(AllRepParticles,1)-1,:) = (AllRepParticles(:,ClosestFrame-2:ClosestFrame+1));
        firstrow = firstrow + size(AllRepParticles,1);
    end
    % convert 0s to nans to show only actively transcribing loci
    AllParticlesArray(AllParticlesArray<=0) = nan;
    meanFluoAcrossThese4Frames = nanmean(AllParticlesArray(:));
    NormAllParticles = AllParticlesArray./meanFluoAcrossThese4Frames;
    
    %Edges = [0:0.05:8];
    Edges = logspace(-1.4,1.1,21);
    %histogram(NormAllParticles,Edges,'Normalization','probability')
    [N,~] = histcounts(NormAllParticles,Edges,'Normalization','probability');
    plot(Edges(1:end-1)+(Edges(2:end)/10),smooth(N),'b-','LineWidth',2);
    set(gca,'xscale','log')
    
end

    
    
end

    







