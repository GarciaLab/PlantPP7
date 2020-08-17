
function SingleLiveExperimentPlots(varargin)

%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

% read arguments
if isempty(varargin)%looks for the folder to analyze
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
    end
end

FilePrefix=[Prefix,'_'];

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
[~, FISHPath, DropboxFolder, MS2CodePath, PreProcPath] = DetermineLocalFolders(Prefix);


%Note that some of this information is redundant given what we get out of
%readMovieDatabase above. We'll have to integrate this better.
%[~,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.csv']);
% ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
% ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
% APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));

% DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
% Dashes=findstr(Prefix,'-');
% PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
%     if isempty(PrefixRow)
%         PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
%         if isempty(PrefixRow)
%             error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
%         end
%     end
%         
% if isempty(PrefixRow)
%     error('Entry not found in MovieDatabase.xlsx')
% end

% ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
% ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};
% APResolution = XLSRaw{PrefixRow,APResolutionColumn};


%Load all the information
HaveCompiled = 0;
if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'CompiledParticles')
    HaveCompiled = 1;
end
if HaveCompiled
    if iscell(CompiledParticles)
        CompiledParticles = CompiledParticles{1};
    end
end

load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
end
if exist([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
end

% % check if a folder with the output of this script exists. If it does,
% % delete it to start fresh.
% ResultsFiguresFolder = ([DropboxFolder,filesep,Prefix,filesep,'ResultsFigures' Prefix]);
% cd DropboxFolder
% if isdir(ResultsFiguresFolder)
%     rmdir(ResultsFiguresFolder)
% end

%Create folder to store the output of this script
ResultsFiguresFolder = ([DropboxFolder,filesep,Prefix,filesep,'ResultsFigures_' Prefix]);
mkdir(ResultsFiguresFolder)

% setup figure parameters
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesFontWeight','bold')

%setup some other stuff
MaxParticleFluo = nanmax([CompiledParticles.Fluo]);
MinParticleFluo = nanmin([CompiledParticles.Fluo]);

MaxMovieTime = (FrameInfo(length(FrameInfo)).Time)/60; %in minutes
MovieTimes = [FrameInfo.Time]/60;

%% Count the number of cells per frame
CellsPerFrame = zeros(1,length(FrameInfo)); 
if exist ('schnitzcells','var')
    for s = 1:length(schnitzcells)
        schnitzFrames = schnitzcells(s).frames;
        CellsPerFrame(schnitzFrames) = CellsPerFrame(schnitzFrames)+1;
    end
    MedianCellNumber = median(CellsPerFrame);
    plot([1:length(FrameInfo)],CellsPerFrame,'-o')
    hold on
    plot([1 length(FrameInfo)],[MedianCellNumber MedianCellNumber],'k-')
    xlabel('frame')
    ylabel('number of cells')
    hold off
    legend ('cells per frame','median')
    saveas(gcf, [ResultsFiguresFolder '\CellsPerFrame.fig'])
else
    MedianCellNumber = length(CompiledParticles);
end


clear schnitzFrames
%% Plot individual traces, one per figure, with and without errorbars
close all
% matrix to store fluo per frame of each particles
AllParticles = nan(max(length(CompiledParticles),round(MedianCellNumber)),length(FrameInfo));
ParticlesPerFrame = zeros(1,length(FrameInfo)); %this counts the number of particles per frame
MinFrames = 0; %minimum number of frames a particle has to have to be considered here

AUPerGFP = 0.076;
GFPPerAU = 1/AUPerGFP;
LoopsPerGFP = 1/2;
PolsPerLoop = 1/48;
PolPerAU = GFPPerAU * LoopsPerGFP * PolsPerLoop;
%AUPerPol = 1/PolPerAU;

%AULabels = [0:100:MaxParticleFluo*1.1];
%NLoopsLabels = AULabels.*LoopsPerAU;
hold on
for p = 1:length(CompiledParticles)
    
        if CompiledParticles(p).Approved > -1
            
        BlankFluo = zeros(1,length(MovieTimes)); % this is to plot a 0 datapoint when the particle is not detected
        BlankError = zeros(1,length(MovieTimes));
        particleFluo = CompiledParticles(p).Fluo;
        particleFluo(isnan(particleFluo)) = 0; 
        particleError = CompiledParticles(p).FluoError;
        particleErrors = ones(1,length(particleFluo)) * particleError;
        particleOffset = CompiledParticles(p).Off;
        particleFrames = CompiledParticles(p).Frame;
        particleTimes = MovieTimes(particleFrames);
        %idx is the index in MovieTimes of the frames where this particle
        %is present
        [~,idxs] = intersect(MovieTimes,particleTimes,'stable'); %indexes of movie frames in which particle is present
        BlankFluo(idxs) = particleFluo; %replace 0s for fluo values when the particle is present
        BlankError(idxs) = particleErrors;   
        
        % add this particle to the arrays where we store the data
        ParticlesPerFrame(particleFrames) = ParticlesPerFrame(particleFrames) + ones(1,length(particleFrames));
        AllParticles(p,:) = BlankFluo;

        if length(particleFrames) > MinFrames
                    
            %plot(MovieTimes,BlankFluo,'ro-','MarkerFaceColor','r','LineWidth',1.5,'MarkerSize',3)
            %hold off
            %shadedErrorBar(MovieTimes,BlankFluo.*LoopsPerAU,BlankError,'lineProps',{'Color',[1 0.5 0.5],'LineWidth',2})
            %yticklabels(num2cell(NLoopsLabels))         
%             NewYticklabels = num2cell(str2double(yticklabels).*LoopsPerAU);
%             yticklabels(NewYticklabels);
            yyaxis left
            errorbar(MovieTimes,BlankFluo,BlankError,'r','CapSize',0,'LineWidth',2)
            ylabel('spot fluorescence (a.u)')
            ylim([0 MaxParticleFluo*1.1])
            yyaxis right
            plot(MovieTimes,BlankFluo.*PolPerAU,'LineStyle','none')
            ylabel('number of transcribing polymerases')
            ylim([ 0 MaxParticleFluo*1.1*PolPerAU])
            
            xlim([0 MaxMovieTime])
            xlabel('time (min)')
            title(['cell #' num2str(p)])
            saveas(gcf, [ResultsFiguresFolder '\' num2str(p) '.fig'])
            close all
        end
    end
end
%save the AllParticles array for the future
save([ResultsFiguresFolder '\AllParticles.mat'],'AllParticles')

clear BlankFluo BlankError particleFluo particleError particleErrors particleOffset particleFrames...
    particleTimes MinFrames
%% Make heatmap of all particles

imagesc(AllParticles);



%% Plot mean transcriptional activity per frame, for just active spots and for all.
close all

% this is for the mean of the actively transcribing ones. 
% in AllParticles.mat *inactive ones are stored as 0*
% we have to convert 0s to nans to get the mean across active particles
% only
AllParticlesWithNans = AllParticles;
AllParticlesWithNans(AllParticlesWithNans==0) = nan;
MeanFluoOn = nanmean(AllParticlesWithNans,1);
SDOn = nanstd(AllParticlesWithNans,1); %standard deviation
SEMOn = SDOn./sqrt(sum(~isnan(AllParticlesWithNans),1)); %standard error of the mean

% now we do it for the mean of active and inactive ones: We create an array called AllParticlesWithOff, which 
% is just AllParticles.mat where we convert nans to zeros so that they count in the calculation.
AllParticlesWithOff = AllParticles;
AllParticlesWithOff(isnan(AllParticlesWithOff)) = 0;
MeanFluoAll = mean(AllParticlesWithOff,1);
SDAll = std(AllParticlesWithOff,1);
SEMAll = SDAll./sqrt(ParticlesPerFrame);

errorbar(MovieTimes,MeanFluoOn,SEMOn,'Color',[0.7 0.9 0.7],...
    'LineStyle','-','LineWidth',1.5,'Marker','o','MarkerFaceColor',[0.4 0.9 0.4],...
    'MarkerEdgeColor','none','CapSize',0)
hold on
errorbar(MovieTimes,MeanFluoAll,SEMAll,'Color',[0.7 0.7 0.9],...
    'LineStyle','-','LineWidth',1.5,'Marker','o','MarkerFaceColor',[0.6 0.6 0.9],...
    'MarkerEdgeColor','none','CapSize',0)
hold off
legend('mean of active loci','mean activity of all loci')
ylabel ('spot fluorescence (a.u)')
xlabel('time (min)')
title(['Mean reporter activity ' Prefix])
ylim([0 400]);
saveas(gcf, [ResultsFiguresFolder '\MeanSpotFluorescence_bars.fig'])
close all

% now the same but with shadedErrorBar
shadedErrorBar(MovieTimes,MeanFluoOn,SDOn./sqrt(ParticlesPerFrame),'lineProps',{'Color',[0.5 0.8 0.5],'LineWidth',2})
hold on
shadedErrorBar(MovieTimes,MeanFluoAll,SDOn./sqrt(ParticlesPerFrame),'lineProps',{'Color',[0.5 0.5 0.8],'LineWidth',2})
hold off

legend('mean of active cells','mean of all cells')
ylabel ('spot fluorescence (a.u)')
xlabel('time (min)')
title(['Mean reporter activity ' Prefix])
ylim([0 400]);
saveas(gcf, [ResultsFiguresFolder '\MeanSpotFluorescence_shaded.fig'])

%% Plot all particles in the same Figure
close all
figure
plot(MovieTimes,AllParticlesWithOff','-','Color',[0.8 0.5 1],'LineWidth',1);
ylabel('time (min)')
xlabel('spot fluorescence')
ylim([0 MaxParticleFluo*1.1])
xlim([0 MaxMovieTime])
title(['All Particles ' Prefix])
saveas(gcf, [ResultsFiguresFolder '\AllParticles.fig'])

%% Plot all particles and the mean of the active ones on top
close all
figure
plot(MovieTimes,AllParticlesWithOff','-','Color',[1 0.7 0.7],'LineWidth',1.2)
hold on
errorbar(MovieTimes,MeanFluoOn,SDOn./sqrt(ParticlesPerFrame),'Color',[0.4 0.9 0.5],...
    'LineStyle','-','LineWidth',1.5,'Marker','o','MarkerFaceColor',[0.1 0.85 0.1],...
    'MarkerEdgeColor','none','CapSize',0)

xlabel('time (min)')
ylabel('spot fluorescence')
ylim([0 max(AllParticles(:))*1.1])
title(['Mean Reporter Activity ' Prefix])
f=get(gca,'Children');
leg = legend([f(end-1),f(1)],'single loci','mean and SEM of active loci');
saveas(gcf, [ResultsFiguresFolder '\AllParticles_MeanOn.fig'])

%% Plot all particles and the mean of all (active and inactive) on top
figure
plot(MovieTimes,AllParticles','-','Color',[1 0.7 0.7],'LineWidth',1.2)
hold on
errorbar(MovieTimes,MeanFluoAll,SDOn./sqrt(ParticlesPerFrame),'Color',[0.4 0.4 0.9],...
    'LineStyle','-','LineWidth',1.5,'Marker','o','MarkerFaceColor',[0.2 0.2 0.9],...
    'MarkerEdgeColor','none','CapSize',0)

xlabel('time (min)')
ylabel('spot fluorescence')
ylim([0 max(AllParticles(:))*1.1])
title(['Reporter activity ' Prefix])
f=get(gca,'Children');
leg = legend([f(end-1),f(1)],'single loci','mean and SEM of all loci');
saveas(gcf, [ResultsFiguresFolder '\AllParticles_MeanAll.fig'])

clear SDOn SEMOn SDAll SEMAll

%% Plot total spot fluorescence
close all
figure
TotalSpotFluoPerFrame = zeros(1,length([FrameInfo.Time]));
for p = 1:length(CompiledParticles)
    %if Particles(CompiledParticles(p).OriginalParticle).Approved == 0;
        ParticleFrames = CompiledParticles(p).Frame;
        ParticleFluo = CompiledParticles(p).Fluo;
        ParticleFluo(isnan(ParticleFluo))=0;
        TotalSpotFluoPerFrame(ParticleFrames) = TotalSpotFluoPerFrame(ParticleFrames)+ParticleFluo;
    %end
end
plot([FrameInfo.Time]./60,TotalSpotFluoPerFrame,'-ob','LineWidth',2,'MarkerFaceColor','b')
xlabel('time (min)')
ylabel('total spot fluorescence')
title(['Total spot fluorescence ' Prefix])
ylim([0 max(TotalSpotFluoPerFrame)*1.2])
saveas(gcf, [ResultsFiguresFolder '\TotalSpotFluorescence.fig'])


%% Number of spots over time
close all
figure
plot([FrameInfo.Time]./60,ParticlesPerFrame,'-ok','LineWidth',2,'MarkerFaceColor','k')
ylabel('number of spots')
xlabel('time (min)')
ylim([0 max(ParticlesPerFrame)*1.2])
title(['Number of spots ' Prefix])
saveas(gcf, [ResultsFiguresFolder '\SpotNumber.fig'])


%% Fraction of active cells over time as: # of cells / # of spots
% ASSUMES 1 NUCLEUS = 1 SPOT
% This is not great since some cells can have more spots than others.
close all
plot([FrameInfo.Time]./60,ParticlesPerFrame./CellsPerFrame,'-ok','LineWidth',2,'MarkerFaceColor','g')
ylabel('number of spots / number of nuclei')
xlabel('time (min)')
ylim([0 1])
title(['# of spots (t) / number of cells (t) ' Prefix])
saveas(gcf, [ResultsFiguresFolder '\FractionON.fig'])
%save([ResultsFiguresFolder '\' num2str(i) 'InstFractionON'],'InstFractionON');

%% Fraction of active cells as # of cells with at least one spot / total # of cells
close all
% first, make an extra field in schnitzcells to track wether or not that
% shnitz has a spot and how many it has. It's a row of as many 0s as frames
% that we'll populate later.
% a second extra field has the identity (row in CompiledParticles) of the
% particle associated with that schnitz.
% we'll make particle-schnitz assignments from scratch based on xy proximity.

for s = 1:length(schnitzcells)
    schnitzFrames = schnitzcells(s).frames;
    schnitzcells(s).NumberOfSpots = zeros(1,length(schnitzFrames));
    schnitzcells(s).ParticleID = zeros(1,length(schnitzFrames));
end

for p = 1:length(CompiledParticles)
    particleFrames = CompiledParticles(p).Frame;
    ParticleXpos = CompiledParticles(p).xPos;
    ParticleYpos = CompiledParticles(p).yPos;
    CompiledParticles(p).SchnitzID = nan(1,length(ParticleFrames));
    for pf = particleFrames
        PartFrameX = ParticleXpos(particleFrames==pf);
        PartFrameY = ParticleYpos(particleFrames==pf);
        Distances = nan(1,length(schnitzcells));
        %find the closest schnitz in this frame to make the assignment
        %frame-wise
        for s = 1:length(schnitzcells)
           schnitzFrames = schnitzcells(s).frames;
           if ismember(pf,schnitzFrames)
               NucFrameX = schnitzcells(s).cenx;
               NucFrameX = NucFrameX(schnitzFrames==pf);
               NucFrameY = schnitzcells(s).ceny;
               NucFrameY = NucFrameY(schnitzFrames==pf);
               Distances(s) = norm(double([NucFrameX NucFrameY]) - [PartFrameX PartFrameY]);
           end
        end
        ClosestNucThisFrame = find(Distances == min(Distances));
        CompiledParticles(p).SchnitzID(particleFrames==pf) = ClosestNucThisFrame;
    end
end

%save CompiledParticles and Schnitzcells

save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_CompiledParts_NucleusID.mat'],'CompiledParticles')


% now that we have nucleus assignments in CompiledParticles, use that info
% to go back to schnitzcells and add it.

for p = 1:length(CompiledParticles)
    ParticleFrames = CompiledParticles(p).Frame;
    ParticleNucleusIDs = CompiledParticles(p).SchnitzID;
    for pfr = ParticleFrames
        NucleusID = (ParticleNucleusIDs(ParticleFrames == pfr));
        NucleusFrames = schnitzcells(NucleusID).frames;
        schnitzcells(NucleusID).NumberOfSpots(NucleusFrames == pfr) = 1 + schnitzcells(NucleusID).NumberOfSpots(NucleusFrames == pfr)
    end
end

% count the number of nuclei with *at least* one spot per frame

NucleiWithSpotPerFrame = zeros(1,length(FrameInfo));
for s = 1:length(schnitzcells)
    NucleusFrames = schnitzcells(s).frames;
    NucleusSpots = schnitzcells(s).NumberOfSpots;
    NucleusSpots(NucleusSpots>1) = 1; %if they have more than one we'll count it as 1 active nucleus
    NucleiWithSpotPerFrame(NucleusFrames) = NucleiWithSpotPerFrame(NucleusFrames) + NucleusSpots;
end

figure
plot([FrameInfo.Time]./60,NucleiWithSpotPerFrame./CellsPerFrame,'k-..')
hold on
plot([FrameInfo.Time]./60,ParticlesPerFrame./CellsPerFrame,'r-.o')
hold off
legend('accounting for multiple spots per nucleus','assuming one nucleus = one particle')
title('Instantaneous fraction on')
xlabel('time (min)')
ylabel('fraction of active cells')
saveas(gcf, [ResultsFiguresFolder '\InstaFractionON.fig'])

InstFractionON = ParticlesPerFrame./CellsPerFrame;

% sanity check: multiplying the mean spot fluorescence of active cells 'MeanON' by
% the number of ative cells 'NucleiWithSpotPerFrame' should be equal to the total spot fluorescence
% 'TotalSpotFluoPerFrame'
figure
plot([FrameInfo.Time]./60,TotalSpotFluoPerFrame,'-bd','MarkerSize',10,'MarkerFaceColor','b')
hold on
plot([FrameInfo.Time]./60,MeanFluoOn .* NucleiWithSpotPerFrame ,'-rs','MarkerFaceColor','r')
hold off
xlabel('time (min)')
ylabel('activity')
legend('total spot fluorescence','mean spot fluorescence of active cells x # of active cells')
saveas(gcf, [ResultsFiguresFolder '\Decomposition_SanityCheck1.fig'])

% Another sanity check: the mean spot fluorescence of all cells
% (active+inactive) 'MeanAll' should be equal to the mean fluorescence of active
% cells 'MeanOn' multiplied by the instantaneous fraction of active cells
% 'InstFractionON'

% ***** Note!!, MeanAll was previosuly defined with a single denominator (the
%( maximum between the number of spots and the median number of cells per
%frame). We should define it with a denominator as a function of time.

MeanFluoAll = TotalSpotFluoPerFrame./CellsPerFrame;
figure
plot([FrameInfo.Time]./60,MeanFluoAll,'-gd','MarkerSize',10,'MarkerFaceColor','g')
hold on
%plot([FrameInfo.Time]./60,TotalSpotFluoPerFrame,'-md','MarkerSize',10,'MarkerFaceColor','m')
plot([FrameInfo.Time]./60,MeanFluoOn .* InstFractionON ,'-ks','MarkerFaceColor','k')
hold off
xlabel('time (min)')
ylabel('activity')
legend('mean spot fluorescence - all cells','mean spot fluorescence of active cells x fraction active cells')
saveas(gcf, [ResultsFiguresFolder '\Decomposition_SanityCheck2.fig'])

save([ResultsFiguresFolder '\InstFractionON'],'InstFractionON');
save([ResultsFiguresFolder '\MeanFluoAll'],'MeanFluoAll');
save([ResultsFiguresFolder '\MeanFluoOn'],'MeanFluoOn');



%clear schnitzFrames particleFrames Xpos Ypos ClosestNucThisFrame Distances NucFrameX NucFrameY PartFrameX PartFrameY
%% Heat map of spot fluorescence and binary map of on/off spots
% AllParticles has nans for undetected frames/spots
% AllParticlesWithOff has zeros instead of nans.
close all
colormap(plasma)
imagesc(AllParticles)
colorbar
title(['All Spots ' Prefix])
saveas(gcf, [ResultsFiguresFolder '\AllParticles_heatmap.fig'])
close all

subplot(1,3,1)
imagesc(~isnan(AllParticles))
colormap(viridis)
title(['All Spots, detected/undetected ' Prefix])

% collapse this figure across columns to show fraction of cells that never
% transcribe vs fraction of cells that transcribed at some point
V = sum(AllParticlesWithOff,2);
V(V>0) = 1; % '1' if the nucleus was ever active, '0' if it was never active
FractionActiveEver = sum(V(V==1))/length(V);
subplot(1,3,2)
imagesc(V)

% collapse across rows to show the instantaneous fraction on
V2 = sum(~isnan(AllParticles),1);
subplot(1,3,3)
imagesc(V2)
colorbar


saveas(gcf, [ResultsFiguresFolder '\AllParticles_ONOFF.fig'])
save([ResultsFiguresFolder '\FractionActiveEver'],'FractionActiveEver');
save([ResultsFiguresFolder '\FractionOn_binary_matrix'],'V');







%% Accumulated mRNA (integrated fluorescence)
% remember that 'AllParticlesWithOff' is an array where off spots have a 0.
% instead, 'AllParticles' has NaNs when a particle is off

TimeVector = [FrameInfo.Time];
%IntegratedFluoAll = trapz(TimeVector,AllParticlesWithOff'); %this counts off spots as 0


for f = 2:length(TimeVector)
    AllSpotsSoFarWithOff = AllParticlesWithOff(:,1:f); %this counts undetected spots as 0
    % remember, the 'total number of potential spots' is the maximum
    % between the median #nuclei(frame) and the number of
    % particles(frame)
    TimeSoFar = TimeVector(1:f)./60;
    IntegralSoFarWithOff(:,f) = trapz(TimeSoFar,AllSpotsSoFarWithOff');
end
subplot(1,2,1)
plot(TimeVector./60,IntegralSoFarWithOff,'-','Color',[.7 .7 .9],'LineWidth',1)
title('Integrated RNA')
xlabel('Time (min)')
ylabel('Integrated spot fluorescence (AU)')

subplot(1,2,2)
histogram(IntegralSoFarWithOff(:,end),'BinMethod','fd','Normalization','probability','FaceColor','b','FaceAlpha',0.5,'EdgeColor','none')
title('Distribution of accumulated RNA')
xlabel('integrated spot fluorescence (AU)')
ylabel('probability')
legend('all cells','detected cells')

saveas(gcf, [ResultsFiguresFolder '\IntegratedFluo_singleCells.fig'])

%% Mean accumulated mRNA per cell, for detected spots and for all spots
close all
SDIntegratedFluoAll = std(IntegralSoFarWithOff,0,1);
SEIntegratedFluoAll = SDIntegratedFluoAll./sqrt(CellsPerFrame);
MeanIntegratedFluoAll = mean(IntegralSoFarWithOff,1);
errorbar([FrameInfo.Time]./60,MeanIntegratedFluoAll,SEIntegratedFluoAll,'Color',[.5 .4 1],...
    'CapSize',0,'LineWidth',2)
title(['Integrated mRNA, all cells' Prefix])
ylabel('mean integrated spot fluorescence')
xlabel('time (min)')
set(gca,'FontSize', 18,'FontWeight','Bold')
saveas(gcf, [ResultsFiguresFolder '\MeanAccumulatedFluoAll.fig'])


%now convert 0s to nans in IntegralSoFarWithOff to integrate over the detected
% cells only
idx = ~IntegralSoFarWithOff;
IntegralSoFarOn = IntegralSoFarWithOff;
IntegralSoFarOn(idx) = nan;
SDIntegratedFluoOn = nanstd(IntegralSoFarOn,1);
SEIntegratedFluoOn = SDIntegratedFluoOn./sqrt(ParticlesPerFrame);
MeanIntegratedFluoOn = nanmean(IntegralSoFarOn,1);
errorbar([FrameInfo.Time]./60,MeanIntegratedFluoOn,SEIntegratedFluoOn,'Color',[1 .5 .4],...
    'CapSize',0,'LineWidth',2)
title(['Integrated mRNA, on cells' Prefix])
ylabel('mean integrated spot fluorescence')
xlabel('time (min)')
set(gca,'FontSize', 18,'FontWeight','Bold')
saveas(gcf, [ResultsFiguresFolder '\MeanAccumulatedFluoOn.fig'])

%save the single loci data arrays
save([ResultsFiguresFolder '\IntegralSoFarWithOff'],'IntegralSoFarWithOff');
save([ResultsFiguresFolder '\IntegralSoFarOn'],'IntegralSoFarOn');


clear SDIntegratedFluoAll SEIntegratedFluoAll MeanIntegratedFluoAll idx IntegralSoFarOn...
    SDIntegratedFluoOn SEIntegratedFluoOn MeanIntegratedFluoOn
%% Mean Offset over time
close all
OffSets = nan(length(CompiledParticles),length(FrameInfo));
for p = 1:length(CompiledParticles)
    PartFrames = CompiledParticles(p).Frame;
    PartOff = CompiledParticles(p).Off;
    OffSets(p,PartFrames) = PartOff;
end

MeanOff = nanmean(OffSets,1);
SDOff = nanstd(OffSets,1);
SEOff = SDOff./ParticlesPerFrame;
errorbar([FrameInfo.Time]./60,MeanOff,SEOff)
title(['Mean Offset ' Prefix])
ylabel('spot offset intensity (a.u)')
xlabel('time (min)')
saveas(gcf, [ResultsFiguresFolder '\MeanOffset.fig'])
clear MeanOff SDOff SEOff PartFrames PartOff OffSets

%% Dynamic range decomposition: mean spot fluorescence and number of spots
% I'll calculate the fold change from the time of first detection. This is
% because some reporters are not detected in the absence of induction so
% the fold change with respoect to t=0 would be infinite.
close all
DetectionThresh = 10;%AU of the weakes possible spot
TotFluoTheo = CellsPerFrame .* MeanFluoOn;
FirstDetectionFrame = min([CompiledParticles.Frame]);
InitialTotFluo = TotalSpotFluoPerFrame(FirstDetectionFrame);
InitialTotFluoTheo = TotFluoTheo(FirstDetectionFrame);
InitialFractionON = InstFractionON(FirstDetectionFrame); 
InitialMeanFluoActive = MeanFluoOn(FirstDetectionFrame);

%normalize to the first detection frame
NormTotFluo = TotalSpotFluoPerFrame./InitialTotFluo;
NormTotFluoTheo = TotFluoTheo./InitialTotFluoTheo;
NormFractionON = InstFractionON./InitialFractionON;
NormMeanFluoActive = MeanFluoOn./InitialMeanFluoActive;

%calculate the fold change with respect to the first detection frame, 
% at the frame of the maximum total fold change
figure
[FoldChangeTotFluo MaxFoldChangeFrame] = nanmax(NormTotFluo);

FoldChangeTotFluoTheo = max(NormTotFluoTheo(FirstDetectionFrame+1:end));%(MaxFoldChangeFrame);
FoldChangeFractionON = max(NormFractionON(FirstDetectionFrame+1:end));%(MaxFoldChangeFrame)
FoldChangeMeanFluoActive = max(NormMeanFluoActive(FirstDetectionFrame+1:end));%(MaxFoldChangeFrame)

B = bar([FoldChangeTotFluo,FoldChangeTotFluoTheo,FoldChangeFractionON,FoldChangeMeanFluoActive],...
    'EdgeColor','none','FaceColor','flat');
hold on
plot([0,5],[1,1],'k-')
hold off
B.CData(1,:) = [0.3 0.3 0.3];
B.CData(2,:) = [0.7 0.7 0.7];
B.CData(3,:) = [1 0.5 0.5];
B.CData(4,:) = [0.5 0.5 1];

ylabel('fold change')
set(gca, 'XTickLabel', {'total spot fluo' ' number x mean' 'number of spots' 'mean spot fluo'})
xtickangle(45)
title(Prefix)
saveas(gcf, [ResultsFiguresFolder '\FoldChanges_bars.fig'])

figure
TimeForPlots = ([TimeVector(FirstDetectionFrame:end)]-TimeVector(FirstDetectionFrame))./60;
plot(TimeForPlots,NormTotFluo(FirstDetectionFrame:end),'k','LineWidth',2)
hold on
plot(TimeForPlots,NormFractionON(FirstDetectionFrame:end),'r','LineWidth',2)
plot(TimeForPlots,NormMeanFluoActive(FirstDetectionFrame:end),'b','LineWidth',2)
plot(TimeForPlots,NormMeanFluoActive(FirstDetectionFrame:end).*NormFractionON(FirstDetectionFrame:end),'g','LineWidth',2)
plot([TimeForPlots(1) TimeForPlots(end)],[1 1],'k-')
hold off
title(Prefix)
ylabel('fold change')
xlabel('time (min)')
legend('total spot fluorescence','fraction of active cells','mean fluorescence of active cells','multiplication of mean x fraction')

saveas(gcf, [ResultsFiguresFolder '\FoldChanges_overtime.fig'])

clear FoldChangeTotFluo MaxFoldChangeFrame FoldChangeTotFluoTheo FoldChangeFractionON...
    FoldChangeMeanFluoActive TotFluoTheo TimeForPlots

%% Distribution of spot fluorescence over time
close all
bins = 20;
edges = linspace(0,ceil(nanmax(AllParticles(:))),bins); %histogram edges
HistData = [];
counter = 1;
PickedFrames = 1:length(FrameInfo)-4; %we'll do a moving average
for frame = PickedFrames;
    [N,edges] = histcounts(AllParticles(:,frame:frame+4),edges);
   % plot(edges(1:end-1),N)
    HistData(:,counter)=N;
    counter = counter+1;
end
% these are the x and y values of a 2d plot, sort of
% the bin edges and the counts in a histogram
x = edges(2:end);%1:bins-1;

%now the y values will correspond to the z values in a 3d plot and the x
%values to the y values
Palette = viridis(size(HistData,2));
counter = 1;
for frame = PickedFrames
    y = smoothdata(HistData(:,frame),3)';
    try
    y = smooth(HistData(:,frame),3)';
    catch
    end
    %Palette(counter,:)
    hFill = fill3(frame*ones(1, bins+1), x([1 1:end end]), [0 y 0],[.9 .9 .9],...
        'LineWidth',1,'FaceAlpha', 1);
    counter = counter+1;
    hold on
end
xlabel('time (min)')
ylabel('spot fluorescence (AU)')
zlabel('number of cells')
title(Prefix)
%FramesToTime = [FrameInfo.Time]./60; %in minutes
%xticks(ceil(FramesToTime(PickedFrames(1:8:end))));
%xticklabels(string(ceil(FramesToTime(PickedFrames(1:8:end)))))
Caz = -90;
Cel = 80;
%v = [-5 -2 5];
view([Caz Cel])
set(gcf, 'Position',  [100, 100, 500, 700])
saveas(gcf, [ResultsFiguresFolder '\AllParticles_joyplot.fig'])
clear bins caz Caz cel Cel counter edges frame HistData i N Palette PickedFrames x X y


%% Label nuclei according to their instantaneous and accumulated spot fluorescence
close all
WekaPath = [PreProcPath, filesep, Prefix, filesep,'tr2dProject\segmentation\weka\'];
ProbabilityMapPath = [WekaPath 'classificationImage0.tif'];
FrameNumber = length(FrameInfo);

% The probability map generated by StartTr2d(Prefix) is stored in
% LivemRNA\Data\PreProcessedData\Prefix\tr2dProject\segmentation\weka
% the file is called
% classificationImage0
% it's a x,y,t tiff file of a maximum projection of the nuclear
% segmentation channel. Intensities represent probability of the pixel
% belonging to the classified class and go from 0 to 1.0.
    
if exist(ProbabilityMapPath) %is there a weka-classified probability map of this movie?
    
    if ~exist([WekaPath 'SmoothNuclearMask.mat']) %was the part that converts the probability map into an array ran before?

        % convert 3D tiff stack (x,y,t) into a 3D matrix
        %originally taken from http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

        FileTif = ProbabilityMapPath;
        InfoImage = imfinfo(FileTif);
        mImage = InfoImage(1).Width;
        nImage = InfoImage(1).Height;
        NumberImages = length(InfoImage);
        FinalImage = zeros(nImage,mImage,NumberImages,'uint16'); %the actual array

        TifLink = Tiff(FileTif, 'r');
        for i=1:NumberImages
           TifLink.setDirectory(i);
           FinalImage(:,:,i)=TifLink.read(); %This is the final 3D mat array of the probability map.
        end
        TifLink.close();

        % Threshold the probability map into a binary mask
        Threshold = 0.5; %above this probability an object corresponds to a nucleus
        NuclearMask = FinalImage>Threshold;

        % look at the images
        % FrameNumber = size(NuclearMask,3);
        % for t = 1:FrameNumber
        %     FrameImage = NuclearMask(:,:,t);
        %     imshow(FrameImage,[])
        %     waitforbuttonpress
        % end

        % filter by feature size
        AllAreas = [];
        threshold = 400; %minimum number of pixels in a nucleus
        for t = 1:FrameNumber
            Image = NuclearMask(:,:,t);
            ImLabel = bwlabel(Image);
            ImFeatures = regionprops(ImLabel);
            for f = 1:length(ImFeatures)
                FeatureArea = ImFeatures(f).Area;
                if FeatureArea < threshold
                    ImLabel(ImLabel==f) = 0;
                end
            end       
            NuclearMask(:,:,t) = ImLabel>0;
        end

        % Fill holes....this doesn't seem to work :/
        % NuclearMask(:,:,t) = imfill(NuclearMask,'holes'); %this takes a few seconds

        % Smooth the edges. Just blur it with conv2() or imfilter(), then threshold again

        windowSize = 6; % Whatever odd integer you want that's more than 1.
        kernel = ones(5)/windowSize^2;
        for t = 1:FrameNumber
            t
            blurredNuclearMask(:,:,t) = conv2(NuclearMask(:,:,t), kernel, 'same');
            SmoothNuclearMask(:,:,t) = blurredNuclearMask(:,:,t) > 0.5;
        end

        % look at the images again to see the results of smoothening
        % for t = 1:FrameNumber
        %     FrameImage = SmoothNuclearMask(:,:,t);
        %     imshow(FrameImage,[])
        %     title(['frame ' num2str(t)])
        %     waitforbuttonpress
        % end

        % clear stuff and save
        save([WekaPath 'SmoothNuclearMask.mat'],'SmoothNuclearMask')
        %clearvars -except SmoothNuclearMask
        
    else %load the previously generated, weka-based nuclear mask
        load([WekaPath 'SmoothNuclearMask.mat']);
    end
           
    % for each frame, go over particles and determine their position
    TotalFrames = length(FrameInfo);

    % Generate the figures

    AllFluo = [CompiledParticles.Fluo];
    MaxParticleFluo = nanmax(AllFluo);
    MinParticleFluo = nanmin(AllFluo);

    myColorMap = plasma((ceil(MaxParticleFluo)*10));
    myColorMap(1,:) = 0.9; %color of the background
    myColorMap(end,:) = 0.7;
    colormap(myColorMap);

    for f = 1:length(FrameInfo)
        f
        BWNuclearMask = bwlabel(SmoothNuclearMask(:,:,f));
        BWNuclearMask = immultiply(BWNuclearMask,0.01); %make values negligible

        for p = 1:length(CompiledParticles)
            p
            ParticleFrames = CompiledParticles(p).Frame;
            ParticleFluo = CompiledParticles(p).Fluo;
            ParticleXPoss = CompiledParticles(p).xPos;
            ParticleYPoss = CompiledParticles(p).yPos;
            if ismember(f,ParticleFrames) 
                ParticleFrameFluo = sum(ParticleFluo(ParticleFrames==f));
                ParticleFramePosX = sum(ParticleXPoss(ParticleFrames==f));
                ParticleFramePosY = sum(ParticleYPoss(ParticleFrames==f));
                NucleusID = BWNuclearMask(ParticleFramePosY,ParticleFramePosX); %if this nucleus already had a particle associated, then we shoudl add up new particles
                if NucleusID == 0 %if the particle is not within a nucleus mask, asign it to the closest
                    [i,j,Values] = find(BWNuclearMask);
                    nonZeroPositions = [i,j];
                    Distances = pdist2(nonZeroPositions,[ParticleFramePosY,ParticleFramePosX]); %euclidean is the default
                    [IndexOfClosest,DistanceToClosest] = find(Distances==min(Distances));
                    NucleusID = Values(IndexOfClosest(1));
                end
                BWNuclearMask(BWNuclearMask==NucleusID) = NucleusID + ParticleFrameFluo;
            end
        end 
        BWNuclearMask(BWNuclearMask<1 & BWNuclearMask>0) = MaxParticleFluo*1.1;
        ColorCodedNuclearMask(:,:,f) = BWNuclearMask;
        caxis([0 MaxParticleFluo])
        %imagesc(BWNuclearMask)%,'AlphaData', .5);
        t = title(['Frame ' num2str(f) ' ' Prefix]);
        set(t, 'Interpreter', 'none')        
        colorbar
        %waitforbuttonpress
    end   
    save([ResultsFiguresFolder '\ColorCodedNuclearMask.mat'],'ColorCodedNuclearMask');
    % now make a figure of the last frame with the accumulated fluorescence
    % this should use interpolation...will deal with it later...
    LastFrameNuclei = SmoothNuclearMask(:,:,end);
    BWLastFrNuclei = bwlabel(LastFrameNuclei);
end


% for p = 1:length(CompiledParticles)
%     
%     ParticleFluo = CompiledParticles(p).fluo;
%     IntegratedParticleFluo = nansum(ParticleFluo);
%     ParticleXPos = CompiledParticles(p).xPos;
%     ParticleYPos = CompiledParticles(p).yPos;
%     XYLastPos = [ParticleXPos(end) ParticleYPos(end)];
%     % now find the closest nucleus in the last frame
%     
% 






% %% CLUSTERING STUFF
% 
% % For each active cell, calculate the average distance to all cells and the
% % average distance to the rest of active cells
% 
% %first, determine which schnitz are present in each frame
% % and also the x and y position of each
% 
% %a matrix that tells whether a nucleus was present (1) or not(0) during
% %each frame
% SchnitzIDsperFrame = zeros(FrameNumber,length(schnitzcells));
% % matrices to store X and Y position of each  nucleus
% SchnitzXPosperFrame = nan(FrameNumber,length(schnitzcells));
% SchnitzYPosperFrame = nan(FrameNumber,length(schnitzcells));
% 
% %populate schintz position matrices
% for fr = 1:FrameNumber
%     for s = 1:length(schnitzcells)
%         if any(schnitzcells(s).frames == fr)
%             if schnitzcells(s).NumberOfSpots(find(schnitzcells(s).frames == fr))
%             SchnitzIDsperFrame(fr,s)= 2;
%             else
%             SchnitzIDsperFrame(fr,s)= 1;
%             end
%             SchnitzXPosperFrame(fr,s) = schnitzcells(s).cenx(find(schnitzcells(s).frames == fr));
%             SchnitzYPosperFrame(fr,s) = schnitzcells(s).ceny(find(schnitzcells(s).frames == fr));           
%         end
%     end
% end
% 
% AvgDistanceToInactive = nan(FrameNumber,length(schnitzcells));
% AvgDistanceToActive = nan(FrameNumber,length(schnitzcells));
% 
% for fr = 1:FrameNumber
%     for s = 1:length(schnitzcells)
% %         if SchnitzIDsperFrame(fr,s)==1 %if this is an inactive cell
% %             SchnitzXPos = SchnitzXPosperFrame(fr,s);
% %             SchnitzYPos = SchnitzYPosperFrame(fr,s);
% %             RestSchnitzXPos = SchnitzXPosperFrame(fr,:);
% %             RestSchnitzYPos = SchnitzXPosperFrame(fr,:);
% %             AvgDistanceToEveryoneElse(fr,s) = nanmean(sqrt((RestSchnitzXPos-SchnitzXPos).^2 + (RestSchnitzYPos-SchnitzYPos).^2));
%          if SchnitzIDsperFrame(fr,s)==2 %if this is an active cell
%             ActiveCellsThisFrame = SchnitzIDsperFrame(fr,:)==2;
%             SchnitzXPos = SchnitzXPosperFrame(fr,s);
%             SchnitzYPos = SchnitzYPosperFrame(fr,s);
%             
%             RestOnSchnitzXPos = SchnitzXPosperFrame(fr,ActiveCellsThisFrame);
%             RestOnSchnitzYPos = SchnitzXPosperFrame(fr,ActiveCellsThisFrame);
%             AvgDistanceToActive(fr,s) = nanmean(sqrt((RestOnSchnitzXPos-SchnitzXPos).^2 + (RestOnSchnitzYPos-SchnitzYPos).^2));            
%             
%             RestOffSchnitzXPos = SchnitzXPosperFrame(fr,~ActiveCellsThisFrame);
%             RestOffSchnitzYPos = SchnitzXPosperFrame(fr,~ActiveCellsThisFrame);
%             AvgDistanceToInactive(fr,s) = nanmean(sqrt((RestOffSchnitzXPos-SchnitzXPos).^2 + (RestOffSchnitzYPos-SchnitzYPos).^2));
%         
%         end
%     end
% end
% 
% % the distance of nucleus to itself is zero. I don't want to take that into account
% AvgDistanceToInactive(AvgDistanceToInactive==0)=nan;
% AvgDistanceToActive(AvgDistanceToActive==0)=nan;
% 
% % plot results
% MicronsPerPixel = FrameInfo(1).PixelSize;
% MeanAverageDistanceAll = MicronsPerPixel*(nanmean(AvgDistanceToInactive'));
% ErrorAverageDistanceAll = nanstd(MicronsPerPixel*(AvgDistanceToInactive'));
% MeanAverageDistanceOn = MicronsPerPixel*(nanmean(AvgDistanceToActive'));
% ErrorAverageDistanceOn = nanstd(MicronsPerPixel*(AvgDistanceToActive'));
% 
% 
% figure
% %errorbar([1:FrameNumber],MeanAverageDistanceAll,ErrorAverageDistanceAll,'b','CapSize',0)
% shadedErrorBar(MovieTimes,MeanAverageDistanceAll,ErrorAverageDistanceAll,'lineProps',{'Color','b','LineWidth',2})
% hold on
% %errorbar([1:FrameNumber],MeanAverageDistanceOn,ErrorAverageDistanceOn,'r','CapSize',0)
% shadedErrorBar(MovieTimes,MeanAverageDistanceOn,ErrorAverageDistanceOn,'lineProps',{'Color','r','LineWidth',2})
% hold off
% ylabel('mean distance (\mum)')
% xlabel('time (min)')
% title({'Average distance all (blue) or active (red)',Prefix})
% saveas(gcf, [ResultsFiguresFolder '\AverageDistance.fig'])
% %close all
% figure
% plot([1:FrameNumber],MeanAverageDistanceAll./MeanAverageDistanceOn ,'k-','LineWidth',3)
% grid on
% ylabel('distance to inactive/distance to active cells')
% xlabel('time (min)')
% title('Ratio of Average Distance of Active Cells')
% saveas(gcf, [ResultsFiguresFolder '\AverageDistanceRatio.fig'])
% %close all
% 
% 
% %% RADIAL FRACTION ACTIVE
% % for each active cell, ask how the fraction of active cells depends on the
% % radial distance to that cell.
% InstantaneousFractionOn_Corrected = (NucleiWithSpotPerFrame-1)./(CellsPerFrame-1);
% % notice I'm correcting the fraction on by subtracting one cell per frame
% % this is to account for the active cell I calculate the radius from
% radii =0:10:sqrt(FrameInfo(1).LinesPerFrame^2 + FrameInfo(1).PixelsPerLine^2);%in pixels! 
% MeanRadialFractionOn = nan(FrameNumber,length(radii));
% SEMRadialFractionOn = nan(FrameNumber,length(radii));
% for fr=1:FrameNumber %find the active nuclei in this frame
%     ThisFrameNuclei = SchnitzIDsperFrame(fr,:);
%     ThisFrameActiveNuclei = ThisFrameNuclei==2; %2 is the label for nuclei associated with a particle
%     NumberOfActiveNuclei = sum(ThisFrameActiveNuclei);
%     ActiveNucleiCounter = 1;
%     FractionActiveMatrix = nan(NumberOfActiveNuclei,length(radii));
%     if NumberOfActiveNuclei > 1
%         for activeNuc = find(ThisFrameActiveNuclei) %loop over active nuclei, calculate distance to all other nuclei
%             ThisNucXPos = SchnitzXPosperFrame(fr,activeNuc);
%             ThisNucYPos = SchnitzYPosperFrame(fr,activeNuc);
%             RestNucXPos = SchnitzXPosperFrame(fr,:);
%             RestNucYPos = SchnitzYPosperFrame(fr,:);
%             DistancesToEveryoneElse = sqrt((RestNucXPos-ThisNucXPos).^2 + (RestNucYPos-ThisNucYPos).^2);            
%             DistancesToEveryoneElse(DistancesToEveryoneElse==0)=nan; % delete distance to self (0)
%             for r = 1:length(radii) %for nuclei within this radius, get the fration active
%                 CloseCells = DistancesToEveryoneElse <= radii(r);
%                 CloseCellsID = SchnitzIDsperFrame(fr,CloseCells);
%                 FractionActiveMatrix(ActiveNucleiCounter,r) = (nansum(CloseCellsID==2)+0)./(length(CloseCellsID)+1);
%                 % that '+1' or '+0+ in the previous line is to count (or not) the reference
%                 % cell in the center, which is by definition active'
%             end
%             ActiveNucleiCounter = ActiveNucleiCounter+1;
%         end
%             FractionActiveMatrix = FractionActiveMatrix/InstantaneousFractionOn_Corrected(fr); %normalize to the overall fraction on in that frame
%             MeanRadialFractionOn(fr,:) = nanmean(FractionActiveMatrix,1);
%             SEMRadialFractionOn(fr,:) = nanstd(FractionActiveMatrix)/sqrt(NumberOfActiveNuclei);
%             clear FractionActiveMatrix
%     end
% end
% 
% 
% figure
% radiiInMicrons = radii*MicronsPerPixel;
% PickedFrames = 1:FrameNumber;
% Palette = viridis(max(PickedFrames));
% for f = PickedFrames
%     %plot(radiiInMicrons,MeanRadialFractionOn(f,:)./InstantaneousFractionOn(f),'-','Color',Palette(f,:),'LineWidth',2)
% 
%     %plot(radiiInMicrons,MeanRadialFractionOn(f,:),'-','Color',Palette(f,:),'LineWidth',2)
% 
%     errorbar(radiiInMicrons,MeanRadialFractionOn(f,:),SEMRadialFractionOn(f,:),...
%         'Color',Palette(f,:),'CapSize',0,'LineWidth',1.5)
% 
% hold on
% end
% plot([radiiInMicrons(1) radiiInMicrons(end)],[1 1],'k-','LineWidth',2)
% hold off
% grid on
% xlabel('distance from active cell (\mum)')
% ylabel('enrichment in fraction active')
% ylim([0 4])
% title({'Radial fraction active',Prefix})
% saveas(gcf, [ResultsFiguresFolder '\RadialFractionActive.fig'])
% 
% 
% %% K-neighbors fraction active
% % for each active cell, ask how the fraction of active cells depends on the
% % radial distance to that cell.
% InstantaneousFractionOn_Corrected = (NucleiWithSpotPerFrame-1)./(CellsPerFrame-1);
% % notice I'm correcting the fraction on by subtracting one cell per frame
% % this is to account for the active cell I calculate the radius from
% MeanNeighborsFractionOn = nan(FrameNumber,max(CellsPerFrame-1));
% SEMNeighborsFractionOn = nan(FrameNumber,max(CellsPerFrame-1));
% for fr=1:FrameNumber %find the active nuclei in this frame
%     NNeighbors =1:CellsPerFrame(fr)-1;
%     ThisFrameNuclei = SchnitzIDsperFrame(fr,:);
%     ThisFrameActiveNuclei = ThisFrameNuclei==2; %2 is the label for nuclei associated with a particle
%     NumberOfActiveNuclei = sum(ThisFrameActiveNuclei);
%     ActiveNucleiCounter = 1;
%     FractionActiveMatrix = nan(NumberOfActiveNuclei,length(NNeighbors));
%     if NumberOfActiveNuclei > 1
%         for activeNuc = find(ThisFrameActiveNuclei) %loop over active nuclei, calculate distance to all other nuclei
%             ThisNucXPos = SchnitzXPosperFrame(fr,activeNuc);
%             ThisNucYPos = SchnitzYPosperFrame(fr,activeNuc);
%             RestNucXPos = SchnitzXPosperFrame(fr,:);
%             RestNucYPos = SchnitzYPosperFrame(fr,:);
%             DistancesToEveryoneElse = sqrt((RestNucXPos-ThisNucXPos).^2 + (RestNucYPos-ThisNucYPos).^2);            
%             DistancesToEveryoneElse(DistancesToEveryoneElse==0)=nan; % delete distance to self (0)
%             for k = 1:length(NNeighbors) %for nuclei within this radius, get the fration active
%                 % find the closest 'k' neighbors and calculate the fraction active
%                 SortedDistances = sort(DistancesToEveryoneElse); %increasing, nans at the end
%                 KCloseCells = ismember(DistancesToEveryoneElse,SortedDistances(1:k));        
%                 %CloseCells = DistancesToEveryoneElse <= radii(r);
%                 KCloseCellsID = SchnitzIDsperFrame(fr,KCloseCells);
%                 FractionActiveMatrix(ActiveNucleiCounter,k) = (nansum(KCloseCellsID==2)+0)./(length(KCloseCellsID)+1);
%                 % that '+1' or '+0+ in the previous line is to count (or not) the reference
%                 % cell in the center, which is by definition active'
%             end
%             ActiveNucleiCounter = ActiveNucleiCounter+1;
%         end
%             FractionActiveMatrix = FractionActiveMatrix/InstantaneousFractionOn_Corrected(fr); %normalize to the overall fraction on in that frame
%             MeanNeighborsFractionOn(fr,[1:length(FractionActiveMatrix)]) = nanmean(FractionActiveMatrix,1);
%             SEMNeighborsFractionOn(fr,[1:length(FractionActiveMatrix)]) = nanstd(FractionActiveMatrix)/sqrt(NumberOfActiveNuclei);
%             clear FractionActiveMatrix
%     end
% end
% 
% 
% figure
% %radiiInMicrons = radii*MicronsPerPixel;
% PickedFrames = 1:FrameNumber;
% Palette = magma(max(PickedFrames));
% for f = PickedFrames
%     %plot(radiiInMicrons,MeanRadialFractionOn(f,:)./InstantaneousFractionOn(f),'-','Color',Palette(f,:),'LineWidth',2)
% 
%     %plot(radiiInMicrons,MeanRadialFractionOn(f,:),'-','Color',Palette(f,:),'LineWidth',2)
% 
%     errorbar([1:max(CellsPerFrame-1)],MeanNeighborsFractionOn(f,:),SEMNeighborsFractionOn(f,:),...
%         'Color',Palette(f,:),'CapSize',0,'LineWidth',1.5)
% 
% hold on
% end
% plot([1 max(CellsPerFrame-1)],[1 1],'k-','LineWidth',2)
% hold off
% grid on
% xlabel('K-neighbors away from active cell')
% ylabel('enrichment in fraction active')
% ylim([0 4])
% title({'K-neighbors fraction active',Prefix})
% saveas(gcf, [ResultsFiguresFolder '\KNeighborsFractionActive.fig'])
    














%% Pairwise particle traces comparisons
% Correlation coefficients

% ZeroBackgroundVector = nan(length(FrameInfo),1);
% NanBackgroundVector = nan(length(FrameInfo),1);
% CorrCutoff = 0.8;
% minFrames = 10;
% 
% %save position and correlation coefficient
% 
% CorrMatrix = zeros(length(CompiledParticles));
% DMatrix = nan(length(CompiledParticles));
% NormMatrix = nan(length(CompiledParticles));
%     
    
% for p1 = 1:length(CompiledParticles)
%     
%     trace1Frames  = CompiledParticles(p1).Frame;
%     trace1Fluo = CompiledParticles(p1).Fluo;
%     LongTrace1 = ZeroBackgroundVector;
%     LongXPos1 = NanBackgroundVector;
%     LongYPos1 = NanBackgroundVector;
%     LongTrace1(trace1Frames) = trace1Fluo;
%     LongXPos1(trace1Frames) = CompiledParticles(p1).xPos;
%     LongYPos1(trace1Frames) = CompiledParticles(p1).yPos;
%     
%     for p2 = 1:length(CompiledParticles)
%         
%         trace2Frames  = CompiledParticles(p2).Frame;
%         trace2Fluo = CompiledParticles(p2).Fluo;
%         LongTrace2 = ZeroBackgroundVector;
%         LongTrace2(trace2Frames) = trace2Fluo;
%         LongXPos2 = NanBackgroundVector;
%         LongYPos2 = NanBackgroundVector;
%         LongTrace2(trace2Frames) = trace2Fluo;
%         LongXPos2(trace2Frames) = CompiledParticles(p2).xPos;
%         LongYPos2(trace2Frames) = CompiledParticles(p2).yPos;
%         
%         ParticleDistance = sqrt((LongXPos1-LongXPos2).^2 + (LongYPos1-LongYPos2).^2);
%         MeanParticleDistance = nanmean(ParticleDistance)*FrameInfo(1).PixelSize;
%         DMatrix(p1,p2) = MeanParticleDistance;
%         
%         R = corrcoef(LongTrace1,LongTrace2,'rows','complete'); %correlation between traces
%         Corr = R(2);
%         CorrMatrix(p1,p2) = Corr;
%         
%         %N = norm(LongTrace1'-LongTrace2'); %euclidean distance
%         %NormMatrix(p1,p2) = N;
%         
%         if p1 == p2 || length(trace1Frames)<minFrames || length(trace2Frames)< minFrames 
%             CorrMatrix(p1,p2)= nan;
%             DMatrix(p1,p2) = nan;
%             %NormMatrix(p1,p2) = nan;
%         end
%         
%         if Corr > CorrCutoff && p1~= p2
%             plot(LongXPos1,LongYPos1,'-r*')
%             hold on
%             plot(LongXPos2,LongYPos2,'-bo')
%             hold off
%             ylim([0 FrameInfo(1).LinesPerFrame])
%             xlim([0 FrameInfo(1).PixelsPerLine])
%             %waitforbuttonpress
%         end
%         
%         subplot(1,2,1)
%         scatter(LongTrace1,LongTrace2,'ro','MarkerFaceColor','r')
%         hold on
%         P = polyfit(LongTrace1,LongTrace2,1); %first output is the corr coef (slope), second is the y intercept
%         Slope = P(1);
%         Yint = P(2);
%         plot(LongTrace1 , Yint + LongTrace1*Slope,'k')
%         title('Correlation')
%         legend(['linear fit, corr =' num2str(Corr)])
%         xlabel(['Trace ' num2str(p1)])
%         ylabel(['Trace ' num2str(p2)])
%         hold off
%         subplot(1,2,2)
%         plot([1:length(LongTrace1)],LongTrace1,'g-*')
%         hold on
%         plot([1:length(LongTrace2)],LongTrace2,'b-o')
%         hold off
%         xlabel('frame')
%         ylabel('spot fluorescence')
%         legend(['trace ' num2str(p1)], ['trace ' num2str(p2)])
%         title(['Transcription traces, distance (\mum)=' num2str(MeanParticleDistance)])
%         %waitforbuttonpress
%         %[p1 p2]
%     end
% end
% figure(1)
% subplot(1,2,1)
% heatmap(DMatrix,'Colormap',viridis)
% subplot(1,2,2)
% heatmap(CorrMatrix,'Colormap',polarmap(jet))

% figure(2)
% subplot(1,2,1)
% plot(DMatrix,CorrMatrix,'ko')
% xlabel('distance (\mum)')
% ylabel('correlation')
% title('Spatial correlation in transcription')
% saveas(gcf, [ResultsFiguresFolder '\SpatialCorrelation.fig'])


% subplot(1,2,2)
% plot(DMatrix,NormMatrix,'bo')
% xlabel('distance (\mum)')
% ylabel('distance in fluo space')
% title('Spatial distance in transcription')



% figure(4)
% palette = viridis(length(FrameInfo));
% for p = 1:length(CompiledParticles)
%     x = CompiledParticles(p).xPos;
%     y = CompiledParticles(p).yPos;
%     Frames = CompiledParticles(p).Frame;
%     for f = 1:length(Frames)
%         if Frames(f) < length(FrameInfo)
%             plot(x,y,'o','MarkerFaceColor',palette(f,:),'MarkerSize',5,'MarkerEdgeColor','none')
%             hold on
%         end
%     end
% end
% title('time on for first heat shock')
% hold off
% saveas(gcf, [ResultsFiguresFolder '\TimeONinSpace.fig'])


%% Histogram of turn on times

% TimeOns = [];
% TimeOnsSecond = [];
% SecondCutOff = 35;
% for p = 1:length(CompiledParticles)
%     [particle,firstFrame] = find(~isnan(AllParticles(p,:)),1);
%     TimeOns(p) = firstFrame;
%     if find(~isnan(AllParticles(p,(SecondCutOff:end))),1) %see if the particle comes up again
%         [particle,firstFrameSecond] = find(~isnan(AllParticles(p,(SecondCutOff:end))),1);
%         TimeOnsSecond(p) = firstFrameSecond;
%     end
% end
% [first,stats] = cdfplot(TimeOns);
% firstData = first.YData
% [second,stats] = cdfplot(TimeOnsSecond);
% secondData = second.YData
% 
% subplot(2,2,1)
% hist(TimeOns,20)
% 
% subplot(2,2,2)
% plot(MovieTimes(1:length(firstData)),firstData,'b','LineWidth',2)
% title('cumulative distribution time of first detection')
% 
% subplot(2,2,3)
% hist(TimeOnsSecond,20)
% h = findobj(gca,'Type','patch');
% h.FaceColor = 'r';
% h.EdgeColor = 'w';
% 
% subplot(2,2,4)
% plot(MovieTimes(SecondCutOff:SecondCutOff+length(secondData)-1),secondData,'r','LineWidth',2)
% title('cumulative distribution time of second detection')
% 
% saveas(gcf, [ResultsFiguresFolder '\TimeOnDistribution.fig'])

%% Comparing spot fluorescence of first vs second pulse
% hold on
% for p = 1:length(CompiledParticles)
%     X = nanmedian(AllParticles(p,(1:35)));
%     Y = nanmedian(AllParticles(p,(36:end)));
%     if sum(~isnan(AllParticles(p,(36:end))))== 0 
%         Y = 0
%     end
%     plot(X,Y,'ro','MarkerFaceColor','r')
% end
% hold off
% xlabel('median fluorescence first pulse')
% ylabel('median fluorescence second pulse')



