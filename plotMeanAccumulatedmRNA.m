function ReplicatesData = plotMeanAccumulatedmRNA(DataStruct,GammaVals,samplingTimes,normalizationTime)

%creates a figure showing the mean and standard error in the accumulated
%mRNA as calculated from integrated spot fluorescence and a simulated
%degradation rate.

% DATASTRUCT is a matlab structure containing the data. Is generated by the
% CombineReplicatesFigures script.

% GAMMAVALS is a 1 x n vector array containing degradation rates in units
% of 1/minute. The function will return as many lines in the figure as
% degradation rates provided

% SAMPLINGTIMES is passed on to MeanAccumulatedmRNA function. Contains
% values of time in minutes for interpolation

% T0 is passed to the MeanAccumulatedmRNA function.

% NORMALIZATION if this is not left empty the values will be normalized
% (i.e divided by) the value at a given time in minutes.
if exist('normalizationTime','var')
    [~,sampligTimeIdx] = find(samplingTimes==normalizationTime);
end
figure % without normalization, show several degradation rates
Palette = viridis(length(GammaVals));
counter = 1;
hold on
for gamma = GammaVals
    ForLegend{counter} = num2str(log(2)/gamma); %degradation rate expressed as half life
    Color = Palette(counter,:);
    Data = MeanAccumulatedmRNA(DataStruct,'MeanAccumulatedFluoOn',samplingTimes,gamma);
    
    if exist('normalizationTime','var')
        Data = Data./nanmean(Data(:,sampligTimeIdx));
    end
    
    ReplicatesData(:,:,counter) = Data; %first dim = replicates, second dim = time (specified by SAMPLINGTIMES argument), third dim = simulated degradation rate
    Ymean = nanmean(Data);
    YError = nanstd(Data)./sqrt(size(Data,1));
    errorbar(samplingTimes,Ymean,YError,'Color',Color,'LineWidth',3,'CapSize',0)
    counter = counter+1;
end
L = legend(ForLegend);
title(L,'half life (min)')
ylabel('accumulated spot fluorescence (mean across replicates)')
xlabel('time (min)')
title('Mean Accumulated Spot Fluorescence')