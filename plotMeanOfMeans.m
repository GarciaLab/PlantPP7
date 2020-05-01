function AllDataArray = plotMeanOfMeans(Struct,Field,PolPerAU)

% this function plots the mean and standard error across replicates.
% replicates correspond to rows in the structure STRUCT
% the data for each replicate has the form of a 1xn array where n is time
% frames and is contained in the field FIELD.

% it also returns a r x t array where r = number of replicates and t = time
% in minutes.

% % for absolute y-axis
% AUPerGFP = 0.076; %this number comes from the absolute calibration
% GFPPerAU = 1/AUPerGFP;
% LoopsPerGFP = 1/2;
% PolsPerLoop = 1/48;
% PolPerAU = GFPPerAU * LoopsPerGFP * PolsPerLoop;

maxMinutes = 200; %unlikely that a movie is longer than this
AllDataArray = nan(length(Struct),maxMinutes); %to store all data

for rep = 1:length(Struct)
    repData = getfield(Struct(rep),Field);
    repAbsTime = Struct(rep).AbsTime;
    
    %linear interpolation to take means across discrete time data
    interpPoints = 500; %number of elements in the interpolated vector
    interpTimeVector = linspace(0,ceil(max(repAbsTime)),interpPoints);
    interpData = interp1(repAbsTime,repData,interpTimeVector);
    
    %fill AllDatarray with the interpolated data
    lastTimePoint = ceil(max(repAbsTime));
    for t = 1:lastTimePoint
        [~,Idx] =find(abs(interpTimeVector-(t-1))==min(abs(interpTimeVector-(t-1))));%find the interpolated timepoint closest to t
        InterpFieldVal = interpData(Idx(1));
        AllDataArray(rep,t) = InterpFieldVal;
    end
end

LastTimePointMinute = find(sum(~isnan(AllDataArray),1),1,'last');
TimeForPlot = [0: LastTimePointMinute-1];
Mean = nanmean(AllDataArray);
Mean = Mean(1:LastTimePointMinute);
SD = nanstd(AllDataArray);
SD = SD(1:LastTimePointMinute);
SE = SD/sqrt(rep);

figure

if contains(Field,'Fluo')
    yyaxis left
    shadedErrorBar(TimeForPlot,Mean,SE,'lineProps',{'Color','b','LineWidth',2})
    %errorbar(TimeForPlot,Mean,SE);
    ylim([0 nanmax(Mean)*1.2])
    ylabel('spot fluorescence')
    
    yyaxis right
    plot(TimeForPlot,Mean*PolPerAU,'LineStyle','none')
    ylabel('elongating RNAP')
    ylim([0 PolPerAU*nanmax(Mean)*1.2])
    
    title(Field)
    xlabel('time since spot detection (min)')
    xlim([0 60])

else
    errorbar(TimeForPlot,Mean,SE);
    shadedErrorBar(TimeForPlot,Mean,SE,'lineProps',{'Color','b','LineWidth',2})
    title(Field)
    xlabel('time since spot detection (min)')
    xlim([0 60])


end

    
