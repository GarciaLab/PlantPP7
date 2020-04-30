function NewStruct = GetMatFiles(DatasetsStruct,Prefixes,DynamicsResultsPath,varargin)

NewStruct = DatasetsStruct;
for i = 1:length(varargin)
    
    for j = 1:length(Prefixes)
        PrefixName = Prefixes{j};
        SingleTracesPath = [DynamicsResultsPath '/' PrefixName '/' ['ResultsFigures_' PrefixName]];
        
        MatFileName = [varargin{i} '.mat'];
        load(fullfile(SingleTracesPath,MatFileName));
        %DatasetsStruct(j).MeanFluoAll = MeanFluoAll;
        NewStruct = setfield(NewStruct,{j},varargin{i},eval(varargin{i}));
    end
    
end
