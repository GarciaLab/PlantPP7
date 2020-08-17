function     [bootstrappedOneSpotError,bootstrappedTwoSpotError,...
    ExpectedTwoSpot,ExpectedErrorTwoSpot,...
    Predicted_p,ErrorPredicted_p] = ...
        MultinomialTestOfTwoSpots_snapshot(FractionOneSpot,FractionTwoSpot,TotalNuclei)

    
% bootstrap the error in the fraction of cells with one spot.
cObsP = @(x) (sum(x)/length(x)); %bootstrapped function: fraction
nSamples = 1000;
OneSpot = round(FractionOneSpot * TotalNuclei);
NoneSpot = [ones(1,OneSpot) zeros(1,TotalNuclei-OneSpot)]; %this is the original sample
[bootObsNOne,~] = bootstrp(nSamples, cObsP,NoneSpot); % this is the bootstrapped sample
bootObsNOne(bootObsNOne>0.5) = nan; %
bootstrappedOneSpotError = nanstd(bootObsNOne);


% bootstrap the error in the number of cells with two spot.
cObsP = @(x) sum(x); %bootstrapped function
nSamples = 1000;
TwoSpot = round(FractionTwoSpot * TotalNuclei);
NoneSpot = [ones(1,TwoSpot) zeros(1,TotalNuclei-TwoSpot)]; %this is the original sample
[bootObsNTwo,~] = bootstrp(nSamples, cObsP,NoneSpot); % this is the bootstrapped sample
bootstrappedTwoSpotError = nanstd(bootObsNTwo/TotalNuclei); %error in the fraction
    
% calculate the predicted p for each bootstrapped observation. That gives a
% distribution of values of p, which has a standard deviation.
bootstrapped_p = 0.5*(1 - sqrt(1-(2.*bootObsNOne)));
Predicted_p = nanmean(bootstrapped_p);
ErrorPredicted_p = nanstd(bootstrapped_p);

bootstrapped_p_sqr = bootstrapped_p.^2;
Var_two_spots = TotalNuclei*bootstrapped_p_sqr.*(1-bootstrapped_p_sqr);
SD_two_spots = sqrt(Var_two_spots);

% % calculate 95% confidence interval in the error estimation
% SEM = std(SD_two_spots)/sqrt(length(SD_two_spots)); % Standard Error
% ts = tinv([0.025  0.975],length(SD_two_spots)-1);% T-Score
% CI = mean(SD_two_spots) + ts*SEM; % Confidence Intervals

% calculate the confidence interval without assuming normality
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
%CIFcn(SD_two_spots,95); %returns lower and upper interval values
ExpectedErrorTwoSpot = CIFcn(SD_two_spots,95);

%ExpectedErrorTwoSpot = max(SD_two_spots)-min(SD_two_spots);

% CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
% ExpectedErrorTwoSpot = CIFcn(bootstrapped_p_sqr,95);

Predicted_TwoSpot = TotalNuclei.*(bootstrapped_p.^2);
ExpectedTwoSpot = nanmean(Predicted_TwoSpot);


%ExpectedErrorTwoSpot = nanstd(Predicted_TwoSpot);

%ExpectedVarianceTwoSpot = TotalNuclei .* (bootstrapped_p.^2) .* (1-(bootstrapped_p.^2));
%ExpectedErrorTwoSpot = nanmean(sqrt(ExpectedVarianceTwoSpot)) + 2*nanstd(sqrt(ExpectedVarianceTwoSpot));


    
% Predicted_p = 0.5*(1 - sqrt(1-(2*FractionOneSpot)));
% ExpectedTwoSpot = TotalNuclei*(Predicted_p^2);
% ExpectedErrorTwoSpot = TotalNuclei * (Predicted_p^2) * (1-(Predicted_p^2)); %this is the variance
% ExpectedErrorTwoSpot = sqrt(ExpectedErrorTwoSpot)*2; % two standard deviations
    
%     
%     
%     
%     
%     
%     
%     
%     
% FractionActivelyTranscribing = FractionOneSpot + FractionTwoSpots;
% 
% 
% 
% %% Binomial distribution - analytical approach
% 
% n = TotalNuclei; %number of independent binomial trials
% valuesForP = [0.05:0.005:0.9];
% counter = 1;
% 
% for p = valuesForP
%     p1 = 2*p*(1-p); %probability of one spot
%     p2 = p^2; %probability of two spots
%     E1 = n * p1; %expected value of cells with one spot
%     E2 = n * p2; %expected value of cells with two spots
% 
%     Var1 = n * p1 * (1-p1); %variance in the number of cells with one spot
%     Var2 = n * p2 * (1-p2); %variance in the number of cells with two spots
%     
%     Error1(counter) = sqrt(Var1);
%     Error2(counter) = sqrt(Var2);
%     ExpectedValues1(counter) = E1;
%     ExpectedValues2(counter) = E2;
% 
%     counter = counter+1;
% end
% 
% FractionOneOrTwo = (ExpectedValues1 + ExpectedValues2)./n;
% %now find the expected number of nuclei with one or two spots
% [~,index] = min(abs(FractionOneOrTwo-FractionActivelyTranscribing));
% 
% ExpectedOneSpot = ExpectedValues1(index);
% ExpectedTwoSpot = ExpectedValues2(index);
% ExpectedErrorOneSpot = 2*Error1(index);
% ExpectedErrorTwoSpot = 2*Error2(index);
% 
% 
% % errorbar(frame,ExpectedOneSpot,ExpectedErrorOneSpot,'r.','CapSize',0)
% % errorbar(frame,ExpectedTwoSpot,ExpectedErrorTwoSpot,'b.','CapSize',0)
% % plot(frame,FractionOneSpot*n,'ro','MarkerFaceColor','r')
% % plot(frame,FractionTwoSpots*n,'bo','MarkerFaceColor','b')
% 
% 
% 
% 
% % 
% % 
% % if bootstrap
% %     % plot data bootstrapping it
% %     [boostrpMeanFrac1Spot,boostrpErrFrac1Spot] = bootstrapFracNSpots(n,ceil(FractionOneSpot*n));
% %     [boostrpMeanFrac2Spot,boostrpErrFrac2Spot] = bootstrapFracNSpots(n,ceil(FractionTwoSpots*n));
% % 
% %     errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac1Spot,'horizontal','k',...
% %         'CapSize',0,'LineWidth',2)
% %     errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac2Spot,'k',...
% %         'CapSize',0,'LineWidth',2)
% %     plot(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,'ko','Color','k','MarkerSize',5,'LineWidth',2)
% %     
% % else
% %     plot(FractionOneSpot,FractionTwoSpots,'ko','Color','k','MarkerSize',10,'LineWidth',2)
% % end
% % 
% % plot(ExpectedValues1./TotalNuclei,ExpectedValues2./TotalNuclei,'r-','LineWidth',2)
% % %hold off
% % 
% % xlabel('nuclei with one spot')
% % ylabel('nuclei with two spots')
% % % ylim([0 1])
% % % xlim([0 1])
% % set(gca,'FontSize',18)
% 
end

