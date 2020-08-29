function [ParticlesMeanFluo,ParticlesDimmestSpotsMean,ParticlesError] =  ...
    FluoErrorDimmestSpots2(Path,Struct,PolPerAU)


ParticlesMeanFluo = []; % the mean fluorescence across frames for each particle
ParticlesDimmestSpotsMean = []; %the mean fluorescence of the three dimmest frames of each particle
ParticlesError =[]; % the fluorescence error of each particle

for rep = 1:length(Struct)
    
    Prefix = Struct(rep).Prefix;
    DynamicsResultsPath = [Path '/' Prefix];
    load([DynamicsResultsPath '/CompiledParticles.mat'],'CompiledParticles');
    if iscell(CompiledParticles) % in old versions of the lab code this could be a cell
        CompiledParticles = CompiledParticles{1};
    end
    
    for p = 1:length(CompiledParticles)
        
        particleFluo = CompiledParticles(p).Fluo;
        particleFluo = particleFluo(particleFluo>1);
        
        nSpots = 3; %number of dimmest spots
        if length(particleFluo)>nSpots
            meanParticleFluo = mean(particleFluo);
            meanParticleDimmestSpots = sort(particleFluo,'ascend');
            ParticleDimmestSpots = meanParticleDimmestSpots(1:nSpots);
            particleError = CompiledParticles(p).FluoError;
            %ParticlesMeanFluo = [ParticlesMeanFluo meanParticleFluo.*PolPerAU];
            ParticlesDimmestSpotsMean = [ParticlesDimmestSpotsMean ParticleDimmestSpots];
            ParticlesError =[ParticlesError particleError.*ones(1,nSpots)];
        end
    end
    
end

% this is to sort particles by fluorescence and extend that sorting to
% the rest of the metrics.
A = [ParticlesDimmestSpotsMean;ParticlesError]';
B = sortrows(A,1,'ascend'); % the second argument specifies which row is the reference one

%% Plot
N = 180;
DimOfDimmestSignal = B((1:N),1);
DimOfDimmestError = B((1:N),2);

DimOfDimmestSignal_RNAP = DimOfDimmestSignal .* PolPerAU;
DimOfDimmestError_RNAP = DimOfDimmestError .* PolPerAU;

close all
figure(2)
hold on
HErr = histogram(DimOfDimmestError_RNAP,'BinWidth',.4,'Normalization','probability','FaceColor','b',...
    'EdgeColor','k')%,'FaceAlpha',.15);
HSign = histogram(DimOfDimmestSignal_RNAP,'BinWidth',.4,'Normalization','probability','FaceColor','g',...
    'EdgeColor','k')%,'FaceAlpha',.15);
xlabel('number of RNAP')
ylabel('frequency')

% % fit to splines
% splineRes = 0.1;
% xE = [HErr.BinEdges(1)-HErr.BinWidth/2 HErr.BinEdges(1:end-1)+HErr.BinWidth/2 HErr.BinEdges(end)];
% yE = [0 HErr.Values 0];
% xxE = xE(1):splineRes:xE(end);
% sE = spline(xE,yE,xxE);
% pE = pchip(xE,yE,xxE)
% plot(xxE,sE,'b','LineWidth',2)
% 
% xS = [HSign.BinEdges(1)-HSign.BinWidth/2 HSign.BinEdges(1:end-1)+HSign.BinWidth/2 HSign.BinEdges(end)];
% yS = [0 HSign.Values 0];
% xxS = xS(1):splineRes:xS(end);
% sS = spline(xS,yS,xxS);
% pS = makima(xS,yS,xxS);
% plot(xxS,sS,'g','LineWidth',2)


%fit errors to polynomials
% yE(xE>5)=0;
% p = polyfit(xE,yE,5)
% yErrorFitted = polyval(p,[xE(1):0.01:xE(end)]);
% plot([xE(1):0.01:xE(end)],yErrorFitted)
% hold on
% %fit signal to polynomials
% p = polyfit(xS,yS,5)
% ySignalFitted = polyval(p,[xS(1):0.01:xS(end)]);
% plot([xS(1):0.01:xS(end)],ySignalFitted)



hold off
ylim([0 max(HErr.Values)*1.2])
xlim([0 9])

%figure(3)


% 
% 
% 
% % the indices of the dimmest quartile based on their mean fluo over time
% IndicesOfDimmestParticles = B(ceil(1:length(CompiledParticles)/Ndimmestparticles),2); 
