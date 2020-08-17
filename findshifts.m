function DataStruct = findshifts(DataStruct)

% I start the induction at a different frame in each experiment.
% to calculate means over time I need to line them up by shifting them in time

Prefixes = {DataStruct.Prefix};
% look at the first prefix to find out what shifts to apply

if strfind(Prefixes{1},'HsfA2')
    frameRate = 1; %minutes
    CellStruct = struct2cell(DataStruct); 
    NumberOfFields = size(CellStruct,1);
    FieldNames = fieldnames(DataStruct);
    NumberOfFields = length(FieldNames);
    disp('working on HsfA2 replicates')
    shifts = [8 0 0 19 0];  
    for f = 1:NumberOfFields
        FieldData = DataStruct.(FieldNames{f});
        if size(FieldData,1)==1 && isnumeric(FieldData)%if the data is a vector array of 1 x frames  
            for p = 1:length(DataStruct) 
                FieldData = DataStruct(p).(FieldNames{f});
                AlignedFieldData = [nan(1,shifts(p)) FieldData];
                    if strcmp(FieldNames(f),'AbsTime') %we have to deal with time a little differently
                            AlignedFieldData = [0:shifts(p)-1 FieldData+shifts(p)]; %to make the first time point 0
                    end
                DataStruct(p).(FieldNames{f}) = AlignedFieldData;
            end
        elseif size(FieldData,1)> 1 %deal with the arrays here
            FieldData = DataStruct(p).(FieldNames{f});
            AlignedFieldData = [ones(size(FieldData,1),shifts(p)).*[1:shifts(p)] FieldData];
            DataStruct(p).(FieldNames{f}) = AlignedFieldData;
        end
    end   

elseif strfind(Prefixes{1},'HSP101')
    %shifts = [7 14 16 16 16 0 16 16 16 8]; 
    shifts = [6 4 1 1 1 10 1 1];
    disp('working on HSP101 replicates')
        
    frameRate = 1; %minutes
    CellStruct = struct2cell(DataStruct); 
    NumberOfFields = size(CellStruct,1);
    FieldNames = fieldnames(DataStruct);
    NumberOfFields = length(FieldNames);

    for f = 1:NumberOfFields
        FieldData = DataStruct.(FieldNames{f});
        if size(FieldData,1)==1 && isnumeric(FieldData)%if the data is a vector array of 1 x frames  
            for p = 1:length(DataStruct)              
                FieldData = DataStruct(p).(FieldNames{f});
                AlignedFieldData = FieldData(shifts(p):end);
                    if strcmp(FieldNames(f),'AbsTime') %we have to deal with time a little differently
                            AlignedFieldData = AlignedFieldData - AlignedFieldData(1); %to make the first time point 0
                    end
                DataStruct(p).(FieldNames{f}) = AlignedFieldData;
            end
        elseif size(FieldData,1)> 1 %deal with the arrays here
            FieldData = DataStruct(p).(FieldNames{f});
            AlignedFieldData = FieldData(:,(shifts(p):end));
            DataStruct(p).(FieldNames{f}) = AlignedFieldData;
        end
    end
end

disp('replicates lined up')
end
    
    
   