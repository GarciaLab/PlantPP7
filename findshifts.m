function DataStruct = findshifts(DataStruct)

% I start the induction at a different frame in each experiment.
% to calculate means over time I need to line them up by shifting them in time

Prefixes = {DataStruct.Prefix};
% look at the first prefix to find out what shifts to apply
if strfind(Prefixes{1},'HSP101')
    shifts = [7 14 16 16 16 0 16 16 16 8]; %for HSP101 datasets
    disp('working on HSP101 replicates')
elseif strfind(Prefixes{1},'HsfA2')
    shifts = [8 0 0 19];  %for HsfA2 datasets
    disp('working on HsfA2 replicates')
else
    disp('Something wrong with the time offsets')
end

frameRate = 1; %minutes
CellStruct = struct2cell(DataStruct); % I need this to get the field names and index them later
NumberOfFields = size(CellStruct,1);
FieldNames = fieldnames(DataStruct);
NumberOfFields = length(FieldNames);

for f = 1:NumberOfFields
    FieldData = DataStruct.(FieldNames{f});
    if size(FieldData,1)==1 && isnumeric(FieldData)%if the data is a vector array of 1 x frames  
%        if ~DataStruct(1).shiftedBefore
%            [DataStruct(:).shiftedBefore] = deal(1); %mark the struct to show this code has been run before
            for p = 1:length(DataStruct)
                %this is a vector of 0s that we're appending in front of the shifted datasets
                TimeShiftVector = [0:frameRate:(shifts(p)-frameRate)];
                if strcmp(FieldNames(f),'AbsTime') %we have to deal with time a little differently
                    DataStruct(p).(FieldNames{f}) = [TimeShiftVector  DataStruct(p).(FieldNames{f}) + shifts(p)];
                else %for the rest of the fields (actual data) we just add 0s at the front
                    DataStruct(p).(FieldNames{f}) = [zeros(1,length(TimeShiftVector)) DataStruct(p).(FieldNames{f})];
                end
            end
%        end
    elseif size(FieldData,1)> 1 %deal with the arrays here
        DataStruct(p).(FieldNames{f}) = [zeros(size(DataStruct(p).(FieldNames{f}),1),length(TimeShiftVector))...
            DataStruct(p).(FieldNames{f})];
        %disp('arrays not finished')
    end
end

disp('replicates lined up')
end
    
    
   