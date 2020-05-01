function plotDynamicRanges(DynRangeFraction,DynRangemeanFOn,DynRangemeanFAll)

% takes 3 arrays each containing the dynamic ranges of each replicte for a
% given metric. These are generated by the plotSingleDynamicRange function.

figure 

hold on

barX = [1 1.5 2];
barY = [nanmean(DynRangemeanFAll) nanmean(DynRangeFraction) nanmean(DynRangemeanFOn)];
bar(barX,barY,0.1,'b')

P1 = plot(1,DynRangemeanFAll,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10);
EB1 = errorbar(1,nanmean(DynRangemeanFAll),std(DynRangemeanFAll)./sqrt(length(DynRangemeanFAll)),...
    'k','CapSize',5,'LineWidth',2);


P2 = plot(1.5,DynRangeFraction,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10);
EB2 = errorbar(1.5,nanmean(DynRangeFraction),std(DynRangeFraction)./sqrt(length(DynRangeFraction)),...
    'k','CapSize',5,'LineWidth',2);


P3 = plot(2,DynRangemeanFOn,'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',10);
EB3 = errorbar(2,nanmean(DynRangemeanFOn),std(DynRangemeanFOn)./sqrt(length(DynRangemeanFOn)),...
    'k','CapSize',5,'LineWidth',2);


plot([0.8 2.2],[1 1],'k-')

hold off
ylim([0 max([DynRangeFraction,DynRangemeanFOn,DynRangemeanFAll])*1.2])
xlim([0.8  2.2])
ylabel('dynamic range')
xticks([1 1.5 2])
xticklabels({'fluo all','fraction on','fluo on'})
xtickangle(45)