function plotFractionCompetent(Struct)

AllData = [];

for p = 1:length(Struct)
    AllParticles = Struct(p).AllParticles;
    %sum across columns to get a vector that will have a number>0 if the
    %particle had fluorescence at some point
    ParticlesCollapsedInTime = sum(AllParticles,2); 
    ParticlesCollapsedInTime(isnan(ParticlesCollapsedInTime)) = 0;
    Competent = ParticlesCollapsedInTime>0;
    FractionCompetent = sum(Competent) / size(Competent,1);
    AllData(p) = FractionCompetent;
end

hold on
P = plot(1,AllData,'ro','MarkerFaceColor','r','MarkerSize',10);
EB = errorbar(1,nanmean(AllData),std(AllData)./sqrt(length(AllData)),...
    'ko','CapSize',0,'MarkerSize',10,'MarkerFaceColor','k',...
    'LineWidth',2);
hold off
ylim([0 1.1])
xlim([0.9 1.1])
ylabel('fraction competent')
legend([P(1) EB(1)],'individual replicates','mean +- SE')
