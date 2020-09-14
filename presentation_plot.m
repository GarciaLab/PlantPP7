CompiledParticles = CompiledParticles{1};
%
clearvars -except CompiledParticles FrameInfo

% %%
% AbsTime = [FrameInfo.Time]./60;
% for frame = 20;%1:40;length(FrameInfo)
%     figure(frame)
% 
%     hold on
%     for p = [1:2 4:length(CompiledParticles)]
%         TimeForPlot = AbsTime(1:frame);
%         FluoForPlot = zeros(1,length(TimeForPlot));
%         ErrorForPlot = zeros(1,length(TimeForPlot));
%     
%         particleFrames = CompiledParticles(p).Frame;
%         particleFluo = CompiledParticles(p).Fluo;
%         particleError = CompiledParticles(p).FluoError;
%         particleError = ones(1,length(particleFluo)).*particleError;
%         
%         upToThisFrame = find(particleFrames>frame,1,'first');
%         upToThisFrame = min([upToThisFrame-1 particleFrames(end)]);
%         [~,FrameIdx] = min(abs(particleFrames-upToThisFrame));
%         %FrameIdx = find(particleFrames==upToThisFrame,1,'first');
%       
%         particlePresentFrames = particleFrames(1:FrameIdx);
%         FluoForPlot(particlePresentFrames) = particleFluo(1:FrameIdx);
%         ErrorForPlot(particlePresentFrames) = particleError(1:FrameIdx);
%         
%         %errorbar(TimeForPlot,FluoForPlot,ErrorForPlot,'r','CapSize',0,'LineWidth',2)
%         plot(AbsTime(particlePresentFrames),FluoForPlot(2:end-1),'r','LineWidth',1.5)
%         %waitforbuttonpress
%     end
%     hold off
%     title([num2str(round(AbsTime(frame))) ' min'])
% end
% 
% 
% %% Plot individual traces, one per figure, with and without errorbars
% close all
% MovieTimes = [FrameInfo.Time]/60;
% % matrix to store fluo per frame of each particles
% %AllParticles = nan(max(length(CompiledParticles),round(MedianCellNumber)),length(FrameInfo));
% %ParticlesPerFrame = zeros(1,length(FrameInfo)); %this counts the number of particles per frame
% MinFrames = 0; %minimum number of frames a particle has to have to be considered here
% 
% AUPerGFP = 0.076;
% GFPPerAU = 1/AUPerGFP;
% LoopsPerGFP = 1/2;
% PolsPerLoop = 1/24;
% PolPerAU = GFPPerAU * LoopsPerGFP * PolsPerLoop;
% 
% 
% 
% for endFrame = [1,10,30,50]
%     MovieTimes = MovieTimes(1:endFrame);
%     hold on
%     for p = 1:length(CompiledParticles)
% 
%             if CompiledParticles(p).Approved > -1
% 
%             BlankFluo = zeros(1,length(MovieTimes)); % this is to plot a 0 datapoint when the particle is not detected
%             BlankError = zeros(1,length(MovieTimes));
%             particleFluo = CompiledParticles(p).Fluo;
%             particleFluo(isnan(particleFluo)) = 0; 
%             particleError = CompiledParticles(p).FluoError;
%             particleErrors = ones(1,length(particleFluo)) * particleError;
%             particleFrames = CompiledParticles(p).Frame;
%             particleTimes = MovieTimes(particleFrames);
%             %idx is the index in MovieTimes of the frames where this particle
%             %is present
%             [~,idxs] = intersect(MovieTimes,particleTimes,'stable'); %indexes of movie frames in which particle is present
%             BlankFluo(idxs) = particleFluo; %replace 0s for fluo values when the particle is present
%             BlankError(idxs) = particleErrors;   
% 
%             % add this particle to the arrays where we store the data
%             ParticlesPerFrame(particleFrames) = ParticlesPerFrame(particleFrames) + ones(1,length(particleFrames));
%             AllParticles(p,:) = BlankFluo;
% 
%             if length(particleFrames) > MinFrames
% 
%                 yyaxis left
%                 errorbar(MovieTimes,BlankFluo,BlankError,'r','CapSize',0,'LineWidth',2)
%                 ylabel('spot fluorescence (a.u)')
%                 ylim([0 MaxParticleFluo*1.1])
%                 yyaxis right
%                 plot(MovieTimes,BlankFluo.*PolPerAU,'LineStyle','none')
%                 ylabel('number of transcribing polymerases')
%                 ylim([ 0 MaxParticleFluo*1.1*PolPerAU])
% 
%                 xlim([0 MaxMovieTime])
%                 xlabel('time (min)')
%                 title(['cell #' num2str(p)])
%                 saveas(gcf, [ResultsFiguresFolder '\' num2str(p) '.fig'])
%                 close all
%             end
%         end
%     end
% end

%%
AbsTime = [FrameInfo.Time]/60;
PolsPerAU = 0.2741;
goodParticles = [1:2 4:size(AllParticles,1)];
Palette = cbrewer('seq', 'Blues', 10);
Palette = Palette(4:end,:);
Palette =  Palette(randperm(end),:);

ChamberTemp2 = ChamberTemp(1:2:end);
ChamberTemp2 = [ChamberTemp2; 24.1;22.4;22.3;22;22];
for f = 1:32%1:length(Frameinfo)
    figure
    yyaxis right
    plot(ChamberTemp2(1+5:f+5),'-','Color',[1 .6 .6],'LineWidth',3)
    ylim([21 42])
    ylabel(['temperature (' char(176) 'C)'])
    hold on
    count=1;
    for p = [1,10,13,15]%1:size(AllParticles,1)]
        Color = Palette(count,:);
        if sum(AllParticles(p,:))>0
            Error = ones(1,f)*CompiledParticles(p).FluoError;       
%             plot(AbsTime(1:f),AllParticles(p,1:f).*PolsPerAU,'Color',Color,...
%                 'LineWidth',2)
            zeroFrames = AllParticles(p,1:f)==0;
            Error(zeroFrames)=0;
            yyaxis left
            
            errorbar([1:f],AllParticles(p,1:f).*PolsPerAU,Error,'-',...
                'Color',Color,'LineWidth',2,'CapSize',0)
            
            ylabel('number of transcribing RNAP')
            ylim([0 250])


        end
        count = count+1;
    end
    
    hold off
    box on
    xlim([0 32])
    title([num2str(f) ' minutes'])
    xlabel('time (min)')
    set(gca,'FontSize',20,'FontName','Arial')
    saveas(gcf,[num2str(round(AbsTime(f))) '_.tif'])
    %waitforbuttonpress
    close all
end


%%
ChamberTime = [30
60
90
120
150
180
210
240
270
300
330
360
390
420
450
480
510
540
570
600
630
660
690
720
750
780
810
840
870
900
930
960
990
1020
1050
1080
1140
1200
1260
1320
1380
1440
1500
1560
1620
1680
1740
1800
1860
1920
1980
2040
2100
2160
2360
2560
2760
2960
3160
3360
3560
3760
3960
4160];
ChamberTime = ChamberTime./60;
ChamberTemp = [25.2
25.2
25.2
25.4
25.6
26
26.4
27
27.6
28.2
28.9
29.7
30.4
31.3
32.1
32.9
33.8
34.7
35.5
36.2
36.9
37.4
37.9
38.2
38.5
38.7
38.9
39.1
39.2
39.3
39.4
39.6
39.7
39.8
39.9
39.9
40.1
40.2
40.4
40.5
40.6
40.7
40.8
40.8
40.9
41
40.4
39.6
39.3
39.1
38.9
38.6
38.4
38.2
36.7
35.2
33.7
32.2
30.6
29.1
27.6
26.3
25.8
25.5];
plot(ChamberTime(1:2:end),ChamberTemp(1:2:end),'o')
