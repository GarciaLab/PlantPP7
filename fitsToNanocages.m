function [FittedMean,FittedSD] = fitsToNanocages(FluorescenceVector,FluoFactor)

% Takes as an argument a 1xn vector of fluorescence values where
% n is the number of observations, i.e the number of nanocages measured.
% The second argument indicates how many times more powerful the laser was
% while imaging the nanocages compared to regular PP7 imaging. Should be 3
% for the 12mers and 5 for the 60mers.

% returns two of the fitted gaussian parameters, the first one is the
% center of the distribution (the mean), the second is the spread of the
% curve (standard deviation).

%figure
H = histogram(FluorescenceVector./FluoFactor,'Normalization','probability',...
    'DisplayStyle','stairs','LineWidth',2,'EdgeColor','b','BinWidth',1.1);
hold on
distributionCenters = mean([H.BinEdges(1:end-1);H.BinEdges(2:end)]);

DataX = distributionCenters; %fluorescence bins
DataY = H.Values; %histogram counts in each bin

% make a gaussian symbolic function to minimize
% this left side of the subtraction is a gaussian evaluated in the domain of the real data
% (DataX), the right hand side is the y values in the data (DataY)
Gaussfun = @(z)z(1) * exp(-(DataX-z(2)).^2./(2*(z(3)^2))) + z(4) - DataY; 

% make guesses about gaussian parameters z(1), z(2), z(3) and z(4) (amplitude, center, spread and baseline)
z0 = [0.02,10,20,0];
lb = [0.01,0,1,0]; % define lower bounds for each parameter
ub = [0.5,20,30,0]; % define upper bounds for each parameter, note that I'm forcing the baseline to be 0

Guess = lsqnonlin(Gaussfun,z0,lb,ub); %fit by minizing the subtraction between actual Y data and a gaussian
FitGauss = Guess(1)*exp(-(DataX-Guess(2)).^2./(2*(Guess(3)^2))) + Guess(4) ;

%function output
FittedMean = Guess(2);
FittedSD = Guess(3);

plot(DataX,FitGauss,'r','LineWidth',2)
%hold off
xlabel('fluorescence (a.u)')
ylabel('frequency')


end
