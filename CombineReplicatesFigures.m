%% Initialize stuff

%get data from figures stored in Single_Traces_Prefix
% make a list of prefixes to use to make paths to the data/figures.

%DynamicsResultsPath = 'S:\Simon\Dropbox\DynamicsResults';
DynamicsResultsPath = '/Users/simon_alamos/Dropbox/DynamicsResults';
CodeRepoPath = '/Users/simon_alamos/Documents/MATLAB/PlantPP7';
RTqPCR_Data_path = '/Users/simon_alamos/Dropbox/PlantPP7_data_for_figures/RTqPCR_raw_data';


% a list of prefixes I want to combine into a single figure

HsfA2Prefixes = {'2019-07-12-13Rb-HsfA2-F-4','2019-04-15-13Rb-HsfA2-F-2','2019-09-14-13Rb-HsfA2-F-5','2019-04-15-13Rb-HsfA2-F'};
EF1alphaPrefixes = {'2019-03-13-EF1alpha-8','2019-04-15-UPG-13Rb-EF1alpha-26','2019-07-12-12R-ef1alpha_17'};

HSP101Prefixes = {'2019-02-16-12R-HSP101_HS_7','2019-02-16-12R-HSP101_HS_5','2019-02-16-12R-HSP101_HS_4',...
    '2019-02-16-12R-HSP101_HS_2','2019-02-16-12R-HSP101_HS_1','2019-02-07-12R-HSP101-HS4',...
    '2019-02-07-12R-HSP101_HS_3','2019-02-07-12R-HSP101_HS_1'};%,'2019-02-07-12R-HSP101_HS2',...
    %'2018-10-23-AL12R-HSP101-3_HS'};

PulsePrefixes = {'2019-02-16-12R-HSP101_HS_RT_HS_4','2019-02-16-12R-HSP101-HS_RT_HS'};

Snapshot2SpotsHSP101Prefixes = {'2020-01-20-AL13Rb-HSP101_9-1-A','2020-01-20-AL13Rb-HSP101_9-1-D',...
    '2020-01-20-AL13Rb-HSP101_9-1-Dc'};

Snapshot2SpotsHsfA2Prefixes = {'2020-01-22-AL13Rb-HsfA2_new-22-1a','2020-01-22-AL13Rb-HsfA2_new-22-1b'};


%for figures
set(0,'defaulttextfontsize',14);
set(0,'defaultaxesfontsize',14);
set(0, 'defaultTextFontName', 'Arial')


%% Run the latest version of the post analysis code here
% GoodPrefixes = [HsfA2Prefixes EF1alphaPrefixes HSP101Prefixes PulsePrefixes];
% GoodPrefixes = HSP101Prefixes;
% 
% 
% tic %to time how long this takes
% for p = 1:length(GoodPrefixes)    
%     PrefixName = GoodPrefixes{p};
%     SingleLiveExperimentPlots(PrefixName);
%     %SpatialClustering(PrefixName);
% end
% toc


%% Pick a group of prefixes and intialize a struct to store everything about each replicate
Prefixes = HSP101Prefixes;
for p = 1:length(Prefixes)
    DatasetsStruct(p).Prefix = Prefixes{p};
end

% Add data from individual replicates to the DatasetsStruct structure
displayStructInfo

DatasetsStruct = GetMatFiles(DatasetsStruct,Prefixes,DynamicsResultsPath,'MeanFluoAll',...
    'MeanFluoOn','InstFractionON','AllParticles','IntegralSoFarWithOff');

FigureNames = {'CellsPerFrame','MeanAccumulatedFluoOn','MeanAccumulatedFluoAll'};
for f = 1:length(FigureNames)
    nameOfFigure = FigureNames{f};
    axis = 1; %sometimes the figure has multiple plots in it, this is the index of the plot. I don't need this
    DatasetsStruct = GetYDataFromFigure(DatasetsStruct,Prefixes,DynamicsResultsPath,nameOfFigure,axis);
end

% Finally, add absolute time information to our struct 
for i = 1:length(Prefixes)
    Prefix = Prefixes{i};
    FrameInfo = load([DynamicsResultsPath '/' Prefix '/FrameInfo.mat']);
    FrameInfo = FrameInfo.FrameInfo;
    DatasetsStruct(i).AbsTime = [FrameInfo.Time]./60; %this is in minutes!
end

% *%*%*%*%*%*%*%*%*%*%*%*% TO DO *%*%*%*%*%*%*%*%*%*%*%*%*%*%

% create a folder to store the analysis figures and .mat files.
%% Plot all experiments in one figure
FieldsToPlot = {'MeanAccumulatedFluoAll','MeanFluoAll','MeanFluoOn','InstFractionON'};
for f = 1:length(FieldsToPlot)
    % we are not really using the output of this function right now, just the figure to look at it
    Out = plotAllPrefixes(DatasetsStruct,FieldsToPlot{f}); 
end
clear Out
%% line up datasets in time and look at results again
alignedDatasetsStruct = DatasetsStruct;
alignedDatasetsStruct = findshifts(alignedDatasetsStruct);

FieldsToPlot = {'MeanAccumulatedFluoAll','MeanFluoAll','MeanFluoOn','InstFractionON'};
for f = 1:length(FieldsToPlot)
    plotAllPrefixes(alignedDatasetsStruct,FieldsToPlot{f});
end


%% Calculate the mean accumulated mRNA based on integrated spot fluorescence
% we will assume no degradation and an arbitrary degradation rate.
% also, with and without normalization to t=60

samplingTimes = [0,5,10,15,30,60]; %in minutes, for interpolation in MeanAccumulatedmRNA function
GammaVals = logspace(-4,-1,6); %degradation rate of mRNA in mRNAs/frame or /min. log(2)/gamma = half life in minutes

%figure without normalization, show several degradation rates
AccmRNA_Data = plotMeanAccumulatedmRNA(alignedDatasetsStruct,GammaVals,samplingTimes);

%figure with normalization, show several degradation rates
FluoAccmRNA_normData = plotMeanAccumulatedmRNA(alignedDatasetsStruct,GammaVals,samplingTimes,60);


%% Compare integrated mRNA from microscopy with RT-qPCR data

% first, take a look at the RT-qPCR results.
deltaCT_results_path = [RTqPCR_Data_path '/qPCR_deltaCTs.mat'];
qPCR_Transgene_norm = plotRTqPCRResults(deltaCT_results_path);

FluoAccmRNA_normData = plotMeanAccumulatedmRNA(alignedDatasetsStruct,GammaVals(1),samplingTimes,60);
[Rsq slope] = qPCR_PP7_Regression(qPCR_Transgene_norm,FluoAccmRNA_normData);

% now we can explore how the regression behaves for different degradation
% rates
RsqResults = []; %to store the R^2 metric
slopeResults = []; %to store the slopes of the linear regression fits
NewGammaVals = logspace(-5,0,20);
counter = 1;
for gamma = NewGammaVals
    FluoAccmRNA_normData = plotMeanAccumulatedmRNA(alignedDatasetsStruct,gamma,samplingTimes,60);
    [Rsq slope] = qPCR_PP7_Regression(qPCR_Transgene_norm,FluoAccmRNA_normData);
    RsqResults(counter) = Rsq;
    slopeResults(counter) = slope;
    counter = counter+1;
    close all
end
figure
hold on
yyaxis left
plot(log(2)./NewGammaVals,RsqResults,'k')
ylabel('R^2')

yyaxis right
plot(log(2)./NewGammaVals,slopeResults,'b')
ylabel('fit slope')

xlabel('reporter transgene mRNA half-life (minutes)')
hold off
set(gca,'xscale','log')

%% Absolute calibration of GFP fluorescence using nanocages

Prefixes60mer = {'2019-06-15-G260G_UBC1-CER-60G_16H_5xpower_2','2019-06-15-G260G_UBC1-CER-60G_16H_5xpower_2_2','2019-06-15-G260G_UBC1-CER-60G_16H_5xpower_2_7',...
    '2019-06-15-G260G_UBC1-CER-60G_16H_5xpower_2_9'};
Prefixes120mer = {'2019-06-08-UBC1-CER-120G_3xPower_series2','2019-06-08-UBC1-CER-120G_3xPower_series3',...
    '2019-06-08-UBC1-CER-120G_3xPower_series4'};

Fluo60merCages = getNanoCagesFluorescence(Prefixes60mer,DynamicsResultsPath);
Fluo120merCages = getNanoCagesFluorescence(Prefixes120mer,DynamicsResultsPath);

AUperGFP = AbsoluteCalibration(Fluo60merCages,Fluo120merCages,5,3);

GFPperPCP = 1; % each PCP is fused to one GFP
PCPperPP7 = 2; %two PCP bind each loop
PP7perPol = 24; %each RNAP has one nascent RNA that has 24 pp7 loops
AUperPol = AUperGFP * GFPperPCP * PCPperPP7 * PP7perPol;
PolPerAU = 1/AUperPol;

%% Detection threshold

clearvars -except alignedDatasetsStruct DatasetsStruct Prefixes DynamicsResultsPath CodeRepoPath PolPerAU
% 
% [Errors,Fluos] = FluoErrorDimmestSpots(DynamicsResultsPath,alignedDatasetsStruct);
% length(Fluos)

[MeanFluos,MeanOfDimmest,Error] = FluoErrorDimmestSpots2(DynamicsResultsPath,alignedDatasetsStruct,PolPerAU);

%% Means across replicates of: Instantaneous fraction on, spot fluo (across all and active only)
% keep around in the workspace only the variable we're going to use
clearvars -except alignedDatasetsStruct DatasetsStruct Prefixes DynamicsResultsPath CodeRepoPath PolPerAU

Fields = {'MeanFluoAll','MeanFluoOn','InstFractionON'};
for f = 1:length(Fields)
    currentField = Fields{f};
    plotMeanOfMeans(alignedDatasetsStruct,currentField,PolPerAU)
end

 %just a sanity check: for each replicate, the multiplication of the instantaneous fraction of
 %active cells at frame t and the mean spot fluorescence across active
 %cells in frame t should be equal to the mean spot fluorescence across
 %all cells in frame t.
plotSanityCheck(alignedDatasetsStruct,'InstFractionON','MeanFluoOn','MeanFluoAll')


%% Fraction competent:
% this is the number of nuclei where a spot was detected in at least one
% frame divided by the total number of nuclei.

plotFractionCompetent(alignedDatasetsStruct)

%% Decomposing the dynamic range in the mean fluorescence across cells into:
% the mean fluorescence across active cells and the instantaneous fraction
% of active cells
% by dynamic range we mean the ratio between F(t=1) and F(t>1) where F(t)
% is going to be the three different metrics of activity
T1 =8; % the first time point for 12R-HSP101-3 replicates
T2 = 39; % the second time point for 12R-HSP101-3 replicates
% T1 = 8; % the first time point for 13Rb-HsfA2-F replicates
% T2 = 15; % the second time point for 13Rb-HsfA2-F replicates

DRFraction = plotSingleDynamicRange(alignedDatasetsStruct,'InstFractionON',T1,T2,PolPerAU);
DRmeanFOn = plotSingleDynamicRange(alignedDatasetsStruct,'MeanFluoOn',T1,T2,PolPerAU);
DRmeanFAll = plotSingleDynamicRange(alignedDatasetsStruct,'MeanFluoAll',T1,T2,PolPerAU);

close all
plotDynamicRanges(DRFraction,DRmeanFOn,DRmeanFAll)

%% Single traces variability and behavior












%% Two spot snapshots








%% Deprecated




















% %% Pool together all the particles from all replicates to generate histograms of 
% % fluorescence over time
% clear AllParticlesPerTime DatasetsStruct
% 
% % get the FrameInfo and AllParticles.mat arrays.
% DynamicsResultsPath = 'E:\Simon\LivemRNA\Data\DynamicsResults';
% CombinedParticlesPerFrame = struct('AllSpots',[]);
% 
% for p = 1:length(Prefixes)
%     TimeOffset = shifts(p); %to align datasets in time   
%     Prefix = Prefixes{p};
%     SingleTracesPath = [DynamicsResultsPath '\' Prefix '\SingleTraces_' Prefix]; 
%     
%     Temp = open([DynamicsResultsPath '/' Prefix '/FrameInfo.mat']);
%     FrameInfo = Temp.FrameInfo;
%     MovieTime = [FrameInfo.Time]./60;
%         
%     Temp = open([SingleTracesPath '/AllParticles.mat']);
%     AllParticles = Temp.AllParticles; %rows are particles, columns are frames
%     AllParticles(AllParticles<0)=nan; %there's a few negative values, probably spurious fits
%     
%     %apply the time offset shift to align datasets by concatenating on the left
%     DatasetsStruct(p).TimeInMin = [0:TimeOffset-frameRate MovieTime+TimeOffset];
%     DatasetsStruct(p).AllParticles = [nan(size(AllParticles,1),TimeOffset) AllParticles];
% 
% 
% end
% clear Temp
% %now combine the time-matched particles into a single array
% LongestMovieTime = ceil(nanmax([DatasetsStruct.TimeInMin]));
% AllParticlesPerTime = nan(1000,LongestMovieTime);% an array to populate with instantaneous particle fluorescence
% for t = 1:LongestMovieTime %loop over time, find the frame closest to this time
%     AllParticles_now = [];
%     for p = 1:length(DatasetsStruct)
%         [diff,idx] = min(abs(DatasetsStruct(p).TimeInMin-t)); %find frame closest to this time t
%         % PROBLEM: for movies shortest than 'longestMovieTime' the closest
%         % frame will be the last one!
%         if diff < 2.5
%             AllParticles_now = [AllParticles_now ; DatasetsStruct(p).AllParticles(:,idx)];
%         end
%     end
%     AllParticlesPerTime((1:length(AllParticles_now)),t) = AllParticles_now;
% end
% plot(AllParticlesPerTime')
% 
% 
% %% MEAN SPOT FLUORESCENCE OVER TIME
% AUPerGFP = 0.076;
% GFPPerAU = 1/AUPerGFP;
% LoopsPerGFP = 1/2;
% PolsPerLoop = 1/48;
% PolPerAU = GFPPerAU * LoopsPerGFP * PolsPerLoop;
% Bins = 15;
% 
% %integrate particles over time. I'm doing a simple sum instead of
% %trapezoidal Riemann sum integration
% AllIntegratedParticles = sum((nansum(AllParticlesPerTime,2))>0);
% 
% MeanAcrossTime = nanmean(AllParticlesPerTime(1:AllIntegratedParticles,:),2);
% AbsMeanOverTime = MeanAcrossTime.*PolPerAU;
% 
% figure
% histogram(log10(MeanAcrossTime),Bins,'Normalization','probability','EdgeColor','none','FaceColor','r')
% title('average activity across time')
% xlabel('log_{10}(mean spot fluorescence(AU))')
% 
% figure
% [HistCounts,HistEdges] = histcounts(log10(MeanAcrossTime),Bins,'Normalization','probability');
% %plot([HistEdges(1) mean([HistEdges(1:end-1);HistEdges(2:end)]) HistEdges(end)],[0 HistCounts 0],'r','LineWidth',2)
% fill([HistEdges(1) mean([HistEdges(1:end-1);HistEdges(2:end)]) HistEdges(end)],[0 HistCounts 0],'r','EdgeColor','none')
% xlabel('log_{10}(mean spot fluorescence(AU))')
% ylabel('frequency')
% title('average activity across time')
% 
% figure
% histogram(log10(AbsMeanOverTime),Bins,'Normalization','probability','EdgeColor','none','FaceColor','b')
% xlabel('log_{10}(mean number of polymerases')
% title('Average number of polymerases across time')
% 
% figure
% [HistCounts,HistEdges] = histcounts(log10(AbsMeanOverTime),Bins,'Normalization','probability');
% %plot([HistEdges(1) mean([HistEdges(1:end-1);HistEdges(2:end)]) HistEdges(end)],[0 HistCounts 0],'b','LineWidth',2)
% fill([HistEdges(1) mean([HistEdges(1:end-1);HistEdges(2:end)]) HistEdges(end)],[0 HistCounts 0],'b','EdgeColor','none')
% xlabel('log_{10}(mean number of polymerases)')
% ylabel('frequency')
% title('Average Polymerases across time')
% %% COMBINE THE ACCUMULATED RNA PER PARTICLE FROM EACH DATASET
% % to show distribution of *accumulated mRNA*
% % the relevant files are called:
% %'IntegralSoFarWithOff.mat' which is an array of particles x frames where
% % the fluorescence was integrated using the trapezoid method over time.
% % Nuclei that never turned on are assumed to have a single particle and
% % this particle is assigned a constant fluorescence of 0.
% 
% % 'IntegralSoFarOn.mat' is the same as 'IntegralSoFarWithOff.mat' except
% % that it only contains real particles, i.e active particles that were
% % detected at least in one frame.
% 
% % Obviously, integration over time depends on how long you're integrating
% % for. What I want to do is take the integrated fluorescence at a set time after 
% % the first particle in the replicate turns on.
% 
% DynamicsResultsPath = 'E:\Simon\LivemRNA\Data\DynamicsResults';
% FileName = '\IntegralSoFarWithOff.mat';
% 
% % find the first frame where there's activity and then the frame closest to
% % N minutes after that first activity frame. Then integrate in between.
% 
% IntegrationTime = 70; %minutes after first spot appears.
% AllTimeAllIntegratedParticles = nan(70,1000); % to combine data from all replicates
% % this matrix is time x particles. At each time, the value of a particle
% % corresponds to the integrated fluorescence up to that timepoint.
% % Because datasets are not perfectly synchronous I find the frame that is
% % the closest to a given time.
% %AllIntegratedParticles =[];
% TimeVector = 1:IntegrationTime;
% 
% for time = TimeVector
%     AllIntegratedParticles =[];
%     for p = 1:length(Prefixes)
%         PrefixName = Prefixes{p};
%         SingleTracesPath = [DynamicsResultsPath '/' PrefixName '/' ['SingleTraces_' PrefixName]];
%         load([SingleTracesPath '/' FileName]) % load 'IntegralSoFarWithOff.mat'
%         FirstFrame = find(nansum(IntegralSoFarWithOff,1)>0,1);
%         load([DynamicsResultsPath '\' PrefixName '\FrameInfo.mat']);
%         TimeInMin = [FrameInfo.Time]./60;
%         IntegrationFrame = find((abs(TimeInMin-time)) == (min(abs(TimeInMin-time))));
%         IntegrationFrame = IntegrationFrame(1);
%         IntegratedFluo = IntegralSoFarWithOff(:,IntegrationFrame);
%         AllIntegratedParticles(end+1:end+length(IntegratedFluo)) = IntegratedFluo;
%     end
%     AllTimeAllIntegratedParticles(time,(1:length(AllIntegratedParticles))) = AllIntegratedParticles;
% end
% 
% AllTimeAllIntegratedParticles(AllTimeAllIntegratedParticles<0) = nan;
% figure
% hold on
% CVs = [];
% NonZeroCVs = [];
% counter = 1;
% Palette = viridis(length(TimeVector));
% for time = TimeVector
%     %figure
% %     histogram(log10(AllTimeAllIntegratedParticles(time,:)),20,'DisplayStyle','stairs','EdgeColor',...
% %         Palette(counter,:),'LineWidth',2)
%     histogram(log10(AllTimeAllIntegratedParticles(time,:)),20,'FaceColor',...
%         Palette(counter,:),'EdgeColor','none','BinWidth',0.2)
%     title([num2str(time) 'min'])
%     ylim([0 14])
%     xlim([0 4])
%     ylabel('number of cells')
%     xlabel('log10(accumulated mRNA) (AU)')
%     CurrentParticles = AllTimeAllIntegratedParticles(time,:);
%     CurrentParticles_on = CurrentParticles(CurrentParticles>0);
%     CVs(counter) = nanstd(CurrentParticles)./nanmean(CurrentParticles);
%     NonZeroCVs(counter) = nanstd(CurrentParticles_on)./nanmean(CurrentParticles_on);
%     counter = counter+1;
% end
% hold off
% 
% figure
% hold on
% plot(TimeVector,CVs*100,'b')
% plot(TimeVector,NonZeroCVs*100,'r')
% 
% % CVs from literature (%)
% Stapel2017 = 295.80; %figure 2C
% Raj2006 = 232.35; %figure 6B
% Ietswaart2017 = 20.38; %figure 2H
% Battich2015 = 200;%figure 1E
% 
% plot([0 60],[Stapel2017 Stapel2017],'k')
% plot([0 60],[Raj2006 Raj2006],'g')
% plot([0 60],[Ietswaart2017 Ietswaart2017],'m')
% plot([0 60],[Battich2015 Battich2015],'y')
% 
% hold off
% ylabel('CV (%)')
% xlabel('Time (min)')
% legend('including 0s','only active cells','Stapel2017','Raj2006','Ietswaart2017','Battich2015')
% 
% Bins = 11;
% figure
% AllIntegratedParticles = sum(~isnan(sum(AllTimeAllIntegratedParticles,1)));
% histogram(log10(AllTimeAllIntegratedParticles(end,1:AllIntegratedParticles)),Bins,...
%     'Normalization','probability','EdgeColor','none','FaceColor','g')
% xlabel('log_{10}(integrated mRNA (AU))')
% 
% figure
% [HistCounts,HistEdges] = histcounts(log10(AllTimeAllIntegratedParticles(end,1:AllIntegratedParticles)),...
%     Bins,'Normalization','probability');
% fill([HistEdges(1) mean([HistEdges(1:end-1);HistEdges(2:end)]) HistEdges(end)],[0 HistCounts 0],'g','EdgeColor','none')
% xlabel('log_{10}(integrated mRNA (AU))')
% ylabel('frequency')
% title('accumulated mRNA distribution')
% 
% %% JOYPLOT OF THE DISTRIBUTION OF INSTANTANEOUS PARTICLE FLUORESCENCE
% % now make a moving window and show the distribution of particle
% % fluorescence
% % get the maximum particle fluo across replicates to define the histogram
% % edges
% 
% %AllParticlesPerTime contains all the particles from all replicates aligned
% %in time
% 
% MaxParticleFluo = nanmax(AllParticlesPerTime(:));
% bins = 20;
% LogEdges = logspace(0,log10(ceil(MaxParticleFluo)*2),bins);
% LinEdges = linspace(0,ceil(MaxParticleFluo),bins);
% HistData = [];
% % do the first detection frame manually
% [N,LogEdges] = histcounts(AllParticlesPerTime(:,1),LogEdges);
% HistData(:,1)=N;
% 
% % then do a moving window of K frames
% K=3;
% counter = 1;
% for frame = K:size(AllParticlesPerTime,2)-K
%     [N,LogEdges] = histcounts(AllParticlesPerTime(:,(frame:frame+K-1)),LogEdges);
%     HistData(:,counter)=N;
%     counter = counter+1;
% end
% % these are the x and y values of a 2d plot, sort of
% % the bin edges and the counts in a histogram
% Bins = LogEdges(2:end);%1:bins-1;
% %x = flip(x);
% 
% %now the y values will correspond to the z values in a 3d plot and the x
% %values to the y values
% figure
% %Palette = viridis(length(1:7:size(AllParticlesPerTime,2)-5));
% PickedFrames = 1:5:90;
% Palette = viridis(length(PickedFrames));
% 
% counter = 1;
% for frame = PickedFrames%size(AllParticlesPerTime,2)-5
%     Counts = smoothdata(HistData(:,frame),3)';
%     %y = HistData(:,frame)';
%     %y = flip(y);
%     %Palette(counter,:)
% %     hFill = fill3(frame*ones(1, bins+1), x([1 1:end end]), [0 y 0],[.8 .8 1],...
% %     'LineWidth',1,'FaceAlpha', 1);
%     %for the logarithmic version:
%     Color = Palette(counter,:);
%     hFill = fill3(frame*ones(1, bins+1), log10(Bins([1 1:end end])), [0 Counts 0],Color,...
%         'LineWidth',1,'FaceAlpha', 1);
%     counter = counter+1;
%     hold on
% end
% xlabel('time (min)')
% ylabel('log_{10}(spot fluorescence (AU))')
% zlabel('number of cells')
% %FramesToTime = [FrameInfo.Time]./60; %in minutes
% %xticks(ceil(FramesToTime(PickedFrames(1:8:end))));
% %xticklabels(string(ceil(FramesToTime(PickedFrames(1:8:end)))))
% Caz = -90;
% Cel = 80;
% view([Caz Cel])
% set(gcf, 'Position',  [100, 100, 500, 700])
% 
% figure
% imagesc(HistData)
% set(gca,'YTick',[1:size(HistData,1)],'YTickLabel',string(round(log10(Bins),2)))
% colorbar;
% colormap viridis
% xlabel('frames')
% ylabel('log_{10}(spot fluorescence (AU))')
% 
% % % now look at histograms as heatmaps
% % figure
% % colormap viridis
% % counter = 1;
% % % PickedFrames = 1:5:size(AllParticlesPerTime,2)-4;
% % PickedFrames = 1:5:72-4;
% 
% % figure
% % AbsTime = [FrameInfo.Time]./60;
% % PickedFramesInMinutes = floor(AbsTime(PickedFrames));
% % for frame = PickedFrames
% %     y = smooth(HistData(:,frame),3)';
% %     SmoothHistData(counter,:) = y;
% %     counter = counter+1;
% % end
% % imagesc(SmoothHistData)
% % xlabel('spot fluorescence bin')
% % ylabel('time (frame)')
% % %set(gca,'XTickLabel',PickedFramesInMinutes)
% % yticklabels(cellstr(num2str(PickedFrames)))
% % xticklabels(cellstr(num2str(log10(edges))))
% % colorbar;
% 
% % yticklabels(cellstr(num2str(PickedFrames)))
% % look at distribution as swarmplots
% 
% % figure
% % %example
% % data = {randn(25, 1), randn(100, 1) + 1, randn(300, 1) - 0.5};
% % plotSpread(data, ...
% %     'xNames', {'25 pts', '100 pts', '300 pts'}, ...
% %     'distributionMarkers', {'o', '+', '.'});
% % 
% % data = {};
% % counter = 1
% % for frame = 1:5:size(AllParticlesPerTime,2)-4
% %    data{counter} = AllParticlesPerTime(:,frame);
% %    counter = counter+1;
% % end
% % plotSpread(data)%,'xNames', {'25 pts', '100 pts', '300 pts'});
% 
% 
% 
% %% Combine FractionActiveEver a.k.a competent fraction
% 
% DynamicsResultsPath = 'E:\Simon\LivemRNA\Data\DynamicsResults';
% FileName = 'FractionActiveEver.mat';
% Prefixes = HSP101Prefixes;
% ArrayForFractionActiveEver = nan(1,length(Prefixes));
% for p = 1:length(Prefixes)
%     %figure(p)
%     PrefixName = Prefixes{p};
%     SingleTracesPath = [DynamicsResultsPath '/' PrefixName '/' ['SingleTraces_' PrefixName]];
%     load([SingleTracesPath '/' FileName])
%     ArrayForFractionActiveEver(p) = FractionActiveEver;
% end
% 
% x = zeros(1,length(Prefixes))';
% y = ArrayForFractionActiveEver';
% %y = randn(10,1);
% %beeswarm(x,y,'up','none',4)
% errorbar(0.05,mean(ArrayForFractionActiveEver),std(ArrayForFractionActiveEver),'or','CapSize',0,'LineWidth',1.5,...
%     'MarkerSize',4)
% hold on
% scatter(rand(10,1)./20,ArrayForFractionActiveEver,'MarkerFaceColor','w')
% %beeswarm(x,y,'up','none',4)
% hold off
% ylim([0 1])
% xlim([-0.2 0.2])
% 
% %% Combine FractionActiveEver a.k.a competent fraction
% 
% DynamicsResultsPath = 'E:\Simon\LivemRNA\Data\DynamicsResults';
% FileName = 'instFractionON.mat';
% Prefixes = HSP101Prefixes;
% ArrayForFractionActiveEver = nan(1,length(Prefixes));
% for p = 1:length(Prefixes)
%     %figure(p)
%     PrefixName = Prefixes{p};
%     SingleTracesPath = [DynamicsResultsPath '/' PrefixName '/' ['SingleTraces_' PrefixName]];
%     load([SingleTracesPath '/' FileName])
%     ArrayForFractionActiveEver(p) = max(InstFractionON);
% end
% 
% x = zeros(1,length(Prefixes))';
% y = ArrayForFractionActiveEver';
% %y = randn(10,1);
% %beeswarm(x,y,'up','none',4)
% errorbar(0.05,mean(ArrayForFractionActiveEver),std(ArrayForFractionActiveEver)./sqrt(length(Prefixes)),'or','CapSize',0,'LineWidth',1.5,...
%     'MarkerSize',4)
% hold on
% scatter(rand(length(Prefixes),1)./20,ArrayForFractionActiveEver,'MarkerFaceColor','w')
% %beeswarm(x,y,'up','none',4)
% hold off
% ylim([0 1])
% xlim([-0.2 0.2])
% 
% 
% %% Plot each cell's rank among it's peers over time
% Prefixes = HSP101Prefixes;
% FileName = '';
% DynamicsResultsPath = 'E:\Simon\LivemRNA\Data\DynamicsResults';
% for p = 5%:length(Prefixes)
%     AllThisParticles = DatasetsStruct(p).AllParticles; 
%     AllThisParticles(isnan(AllThisParticles))=0;
%     ThisParticlesBinned = nan(size(AllThisParticles));
%     % the first column is the first frame at which a particle was detected.
%     % rows are different particles
%     ThisTime = DatasetsStruct(p).FrameInfoTime;
%     for frame = 1:size(AllThisParticles,2)
%         ThisFrameParticles = AllThisParticles(:,frame);
%         Edges = 0:(ceil(nanmax(ThisFrameParticles))+1)/20:ceil(nanmax(ThisFrameParticles))+1;
%         Edges = logspace(1,log10(ceil(nanmax(ThisFrameParticles))+1),10);
%         Edges = logspace(1,log10(ceil(nanmax(ThisFrameParticles))+1),35);
%         BinnedFrameParticles = discretize(ThisFrameParticles,Edges);
%         ThisParticlesBinned(:,frame) = BinnedFrameParticles;
%     end
%     
%     figure(p)
%     plot(ThisParticlesBinned','-','LineWidth',1.2)         
%     %ylim([0 21])
%     ylabel('activity bin')
%     xlabel('time (frame)')
%     title(Prefix)
%     set(0, 'DefaultFigureColormap', jet(64))
%     
% end
% 
% 
% Palette = flip(viridis(size(ThisParticlesBinned,1)));
% hold on
% for spot = 1:size(ThisParticlesBinned,1)
%     figure
%     randomNumber = randi([1 size(ThisParticlesBinned,1)])
%     LineColor = Palette(randomNumber,:);
%     plot(ThisTime,ThisParticlesBinned(spot,:),'-','LineWidth',1.2,'Color',LineColor) 
%     hold off
% ylim([0 35])
% ylabel('activity bin')
% xlabel('time (frame)')
% title(Prefix)
% end
% 
% 
%  
% 
% 
% 
% 
% 
% 
