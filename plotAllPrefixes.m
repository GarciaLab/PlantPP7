function DataArray = plotAllPrefixes(Struct,FieldY)

% returns a r x t 2 dimensional array DATAATTAY where rows (r) are
% replicates (sepcified by the number of rows in STRUCT) and columns are
% timepoints.

figure 
hold on
LegendNames = {};
ColorMap = viridis(length(Struct));
maxFrames = 200; %movies won't be longer than this
DataArray = nan(length(Struct),maxFrames);
for i = 1:length(Struct)
    X = Struct(i).AbsTime;
    Y = getfield(Struct(i),FieldY);
    legendNames{i} = ['dataset #' num2str(i) ' ' Struct(i).Prefix];
    plot(X,Y,'LineWidth',3,'Color',ColorMap(i,:))
    DataArray(i,1:length(X)) = Y;
end
hold off
legend(legendNames,'Interpreter', 'none')
xlabel('time (min)')
ylabel('mean accumulated spot fluorescence per cell')
title(FieldY)
end

