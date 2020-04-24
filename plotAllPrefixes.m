function plotAllPrefixes(Struct,FieldX,FieldY)

figure 
hold on
LegendNames = {};
ColorMap = viridis(length(Struct));
for i = 1:length(Struct)
    X = getfield(Struct(i),FieldX);
    Y = getfield(Struct(i),FieldY);
    legendNames{i} = ['dataset #' num2str(i) ' ' Struct(i).Prefix];
    plot(X,Y,'LineWidth',3,'Color',ColorMap(i,:))
end
hold off
legend(legendNames,'Interpreter', 'none')
xlabel('time (min)')
ylabel('mean accumulated spot fluorescence per cell')
title(FieldY)

end

