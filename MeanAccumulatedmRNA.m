function AccFluoDataForMean_deg =  ...
    MeanAccumulatedmRNA(Struct,AccumulatedFluoField,samplingTimes,gamma)
%returs an array 'AccFluoDataForMean_deg' containing the accumulated fluorescence
%(columns) of all replicates (rows) with a simulated degradation rate given
%by gamma


% we are going to take the mean of discrete time-dependent data
% however, the datapoints are not taken at the exact same time
% this means we have to do some sort of interpolation.
% samplingTimes = specific timepoints that we want to compare with RT-qPCR
% gamma = simulated degradation rate of the reporter mRNA in units of mRNAs/min, 
% FYI mRNA half life =  log(2)/gamma 
% Some degradation rates and their corresponding half-lifes:
% gamma = 0.0058; % half life of two hours
% gamma = 0.0023; %half life of five hours
% gamma = 0.0462; % half life of 15 minutes
% gamma = 0.0231; % half life of 30 minutes
% gamma = 0.1386; %half life of 5 minutes

% simulated accumulated mRNA from all replicates will be stored in these vectors
% AccFluoDataForMean_deg(length(Struct),length(samplingTimes)) = [];
% AccFluoDataForMean(length(Struct),length(samplingTimes)) = [];


for p = 1:length(Struct)
    AccuFluoDeg = [] ; %here we will store the accumulated mRNA values adjusted by degradation from our simulation
    Time = Struct(p).AbsTime;
    AccuFluo = getfield(Struct(p),AccumulatedFluoField);
    AccuFluoDeg(1) = AccuFluo(1); %initialize the accumulated mRNA simulation accounting for degradation
    
    %now do the rest of the simulation timepoints in a loop
    for t = 2:length(Time)
        mRNApreviousStep = AccuFluoDeg(t-1); %number of mRNAs from the previous time step
        Production = AccuFluo(t) - AccuFluo(t-1); %mRNA produced in a time step of delta t
        Degradation = -gamma * AccuFluoDeg(t-1); %mRNA degraded in a time step of delta t
        % actuallize the vector where we store the simulation results
        AccuFluoDeg = [AccuFluoDeg nansum([mRNApreviousStep;Degradation;Production])]; %append the simulation results to a growing vector
    end

    % to get the mean across the simulated accumulated mRNA we will first
    % interpolate
    interpPoints = 500; %number of elements in the interpolated vector
    interpVector = linspace(0,ceil(max(Time)),interpPoints);
    %interpAccFluo = interp1(Time,AccuFluo,interpVector);
    interpAccFluoDeg = interp1(Time,AccuFluoDeg,interpVector);
    
    for t = 1:length(samplingTimes)
        [dummy index ] = min(abs(interpVector-samplingTimes(t))); %index of the point closest to the picked sampling time
        %AccFluoDataForMean(p,t) = interpAccFluo(index);
        AccFluoDataForMean_deg(p,t) = interpAccFluoDeg(index);
    end
    
end




% 
% 
% 
