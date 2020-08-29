function [AllParticlesPerTime AllInstaFractionOn] = makeAllParticlesPerTime(alignedDatasetsStruct,baseline)


% 'AllParticlesPerTime' array contains all the particles from all replicates aligned
%in time

%find the longest movie to define the last time point
Durations = [];
for rep = 1:length(alignedDatasetsStruct)
    repAbsTime = alignedDatasetsStruct(rep).AbsTime;
    Durations = [Durations repAbsTime(end)];
end
LastTimePoint = ceil(max(Durations));

% create the otuput array of particles x time that we'll then populate
AllParticlesPerTime = nan(300,LastTimePoint);

    


% align particles in time according to the time of first detection.
% frameRateInMinutes = 1;
% for rep = 1:length(alignedDatasetsStruct)      
%     thisRepAbsTime  = alignedDatasetsStruct(rep).AbsTime;
%     %find the first column with particles
%     thisRepParticles  = alignedDatasetsStruct(rep).AllParticles;
%     thisRepParticles(isnan(thisRepParticles)) = baseline;
%     % this frame is the first frame with particles
%     firstColumnWithParticles = find(sum(thisRepParticles,1)>baseline,1,'first');
%     % the absolute time of this frame will be the new time zero
%     newT0 = thisRepAbsTime(firstColumnWithParticles);
%     alignedDatasetsStruct(rep).AbsTime = alignedDatasetsStruct(rep).AbsTime - newT0 + frameRateInMinutes;
% end

Prefixes = {alignedDatasetsStruct.Prefix};
% look at the first prefix to find out what shifts to apply

if strfind(Prefixes{1},'HSP101')

    for t =1:LastTimePoint

        for rep = 1:length(alignedDatasetsStruct)      
            % get the right particles of this rep in this time point
            [distance,nearestFrame] = min(abs(alignedDatasetsStruct(rep).AbsTime-t)); 

            thisRepParticlesThisFrame = alignedDatasetsStruct(rep).AllParticles;
            thisRepParticlesThisFrame = thisRepParticlesThisFrame(:,nearestFrame);
            thisRepParticlesThisFrame(isnan(thisRepParticlesThisFrame)) = baseline;
            % add them to the output array

            %find the first row that is still a nan in this column
            firstNanRow = find(isnan(AllParticlesPerTime(:,t)),1,'first');           
            AllParticlesPerTime(firstNanRow:firstNanRow+size(thisRepParticlesThisFrame,1)-1,t) = thisRepParticlesThisFrame;       

            %find the instantaneous fraction on this frame


        end
    end

    LastRowWithParticles = find(~isnan(sum(AllParticlesPerTime,2)),1,'last');
    AllParticlesPerTime = AllParticlesPerTime(1:LastRowWithParticles,:);

elseif strfind(Prefixes{1},'HsfA2')
    shifts = [8 0 0 19 0]
    
    Durations = [];
    for rep = 1:length(alignedDatasetsStruct)
        repAbsTime = alignedDatasetsStruct(rep).AbsTime;
        Durations = [Durations repAbsTime(end)];
    end
    LastTimePoint = floor(min(Durations));


    for t =1:LastTimePoint

        for rep = 1:length(alignedDatasetsStruct) 
            thisRepShift = shifts(rep);
        % % get the right particles of this rep in this time point
        [distance,nearestFrame] = min(abs(alignedDatasetsStruct(rep).AbsTime-t));
        
            if distance < 2

                thisRepParticlesThisFrame = alignedDatasetsStruct(rep).AllParticles;
                thisRepParticlesThisFrame = thisRepParticlesThisFrame(:,nearestFrame);
                %thisRepParticlesThisFrame = thisRepParticlesThisFrame(:,nearestFrame);
                %thisRepParticlesThisFrame(isnan(thisRepParticlesThisFrame)) = baseline;
                % add them to the output array

                % %find the first row that is still a nan in this column
                 firstNanRow = find(isnan(AllParticlesPerTime(:,t)),1,'first');           
        %         AllParticlesPerTime(firstNanRow:firstNanRow+size(thisRepParticlesThisFrame,1)-1,t) = thisRepParticlesThisFrame;       

                AllParticlesPerTime(firstNanRow:firstNanRow+size(thisRepParticlesThisFrame,1)-1,t) = thisRepParticlesThisFrame;       

            end

        end
    end

    LastRowWithParticles = find(nansum(AllParticlesPerTime,2)>0,1,'last');
    AllParticlesPerTime = AllParticlesPerTime(1:LastRowWithParticles,:);
    
end


end
