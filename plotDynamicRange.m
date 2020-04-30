function DynamicRanges = plotDynamicRange(Struct,Field,t1,t2)

% returns the ratio of the value in the field FIELD of the struct STRUCT at
% timepoints t2 vs t1.
% output is a 1 x r array called DYNAMICRANGES where r is the
% number of replicates.

AllDataArray = plotMeanOfMeans(Struct,Field);
%close all % I think this function generates a figure


DynamicRanges = nan(1,length(Struct));

for rep = 1:length(Struct)
    repAbsTime = Struct(rep).AbsTime;
    repYData = smooth(getfield(Struct(rep),Field));

    % interpolate time and y data, then find the value of data at the
    % timepoint closest to t1 and t2
    interpPoints = 500; %number of elements in the interpolated vector
    interpTimeVector = linspace(0,ceil(max(repAbsTime)),interpPoints);
    interpData = interp1(repAbsTime,repYData,interpTimeVector);
    
    [~,T1Idx] =find(abs(interpTimeVector-t1)==min(abs(interpTimeVector-t1)));%find the interpolated timepoint closest to t1
    [~,T2Idx] =find(abs(interpTimeVector-t2)==min(abs(interpTimeVector-t2)));%find the interpolated timepoint closest to t1
    
    Ratio = interpData(T2Idx(1))/interpData(T1Idx(1)); %sometimes find can return two indices, that's why I'm indexing the index!
    DynamicRanges(rep) = Ratio;
    
%     figure
%     hold on
%     plot(interpTimeVector,interpData,'b-','LineWidth',2)
%     plot([t1 t1],[0 max(interpData)],'g.-')
%     plot([t2 t2],[0 max(interpData)],'r.-')
%     plot([0 max(interpTimeVector)],[interpData(T2Idx(1)) interpData(T2Idx(1))],'r.-')
%     plot([0 max(interpTimeVector)],[interpData(T1Idx(1)) interpData(T1Idx(1))],'g.-')
%     hold off
%     title(['Prefix' num2str(rep) '; Ratio = ' num2str(Ratio)])
%     
%     waitforbuttonpress

    
end

% this is to see what happens if we first take the mean and the error and
% then calculate the dynamic range. To do this we use the output of 
% plotMeanOfMeans(Struct,Field)
meanYData = nanmean(AllDataArray,1);
seYData = std(AllDataArray,1)./sqrt(length(DynamicRanges));
% now take the ratio propagating the error from each timepoint
propRatio = meanYData(t2)/meanYData(t1);
ErrorT1 = seYData(t1);
ErrorT2 = seYData(t2);
RatioPropError = sqrt((ErrorT2/meanYData(t2))^2 + (ErrorT1/meanYData(t1))^2);




figure 

hold on
P = plot(1,DynamicRanges,'ro','MarkerFaceColor','r','MarkerSize',10);
EB = errorbar(1,nanmean(DynamicRanges),std(DynamicRanges)./sqrt(length(DynamicRanges)),...
    'ko','CapSize',0,'MarkerSize',10,'MarkerFaceColor','k',...
    'LineWidth',2);
PropEB = errorbar(1,propRatio,RatioPropError,'bo','CapSize',0,'MarkerSize',...
    10,'MarkerFaceColor','b','LineWidth',2);
labelpoints(1,DynamicRanges,string([1:length(DynamicRanges)]),'NW',1,0.6)

hold off
ylim([0 max(DynamicRanges)*1.2])
ylabel('dynamic range')
legend([P(1) EB(1) PropEB(1)],'individual replicates','mean +- SE','propagated')
title([Field ' ' num2str(t2) '/' num2str(t1)])