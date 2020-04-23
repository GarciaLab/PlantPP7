function updatedStruct = GetYDataFromFigures(DatasetsStruct,Prefixes,DynamicsResultsPath,figureName,axis)

updatedStruct = DatasetsStruct;
    
for j = 1:length(Prefixes)
    PrefixName = Prefixes{j};
    SingleTracesPath = [DynamicsResultsPath '/' PrefixName '/' ['ResultsFigures_' PrefixName]];

    FigFileName = [figureName '.fig'];
    F = open([SingleTracesPath '/' FigFileName]);
    YdataObjs = findobj(F,'-property','YData');
    FigureYData = YdataObjs(axis).YData;
    updatedStruct = setfield(updatedStruct,{j},figureName,FigureYData);

    close all
end


end

