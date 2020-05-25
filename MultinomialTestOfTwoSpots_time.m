function     [ExpectedOneSpot,ExpectedTwoSpot,...
        ExpectedErrorOneSpot,ExpectedErrorTwoSpot] = ...
        MultinomialTestOfTwoSpots_time(FractionOneSpot,FractionTwoSpots,TotalNuclei,frame)


FractionActivelyTranscribing = FractionOneSpot + FractionTwoSpots;



%% Binomial distribution - analytical approach

n = TotalNuclei; %number of independent binomial trials
valuesForP = [0.05:0.005:0.9];
counter = 1;

for p = valuesForP
    p1 = 2*p*(1-p); %probability of one spot
    p2 = p^2; %probability of two spots
    E1 = n * p1; %expected value of cells with one spot
    E2 = n * p2; %expected value of cells with two spots

    Var1 = n * p1 * (1-p1); %variance in the number of cells with one spot
    Var2 = n * p2 * (1-p2); %variance in the number of cells with two spots
    
    Error1(counter) = sqrt(Var1);
    Error2(counter) = sqrt(Var2);
    ExpectedValues1(counter) = E1;
    ExpectedValues2(counter) = E2;

    counter = counter+1;
end

FractionOneOrTwo = (ExpectedValues1 + ExpectedValues2)./n;
%now find the expected number of nuclei with one or two spots
[~,index] = min(abs(FractionOneOrTwo-FractionActivelyTranscribing));

ExpectedOneSpot = ExpectedValues1(index);
ExpectedTwoSpot = ExpectedValues2(index);
ExpectedErrorOneSpot = 2*Error1(index);
ExpectedErrorTwoSpot = 2*Error2(index);


% errorbar(frame,ExpectedOneSpot,ExpectedErrorOneSpot,'r.','CapSize',0)
% errorbar(frame,ExpectedTwoSpot,ExpectedErrorTwoSpot,'b.','CapSize',0)
% plot(frame,FractionOneSpot*n,'ro','MarkerFaceColor','r')
% plot(frame,FractionTwoSpots*n,'bo','MarkerFaceColor','b')




% 
% 
% if bootstrap
%     % plot data bootstrapping it
%     [boostrpMeanFrac1Spot,boostrpErrFrac1Spot] = bootstrapFracNSpots(n,ceil(FractionOneSpot*n));
%     [boostrpMeanFrac2Spot,boostrpErrFrac2Spot] = bootstrapFracNSpots(n,ceil(FractionTwoSpots*n));
% 
%     errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac1Spot,'horizontal','k',...
%         'CapSize',0,'LineWidth',2)
%     errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac2Spot,'k',...
%         'CapSize',0,'LineWidth',2)
%     plot(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,'ko','Color','k','MarkerSize',5,'LineWidth',2)
%     
% else
%     plot(FractionOneSpot,FractionTwoSpots,'ko','Color','k','MarkerSize',10,'LineWidth',2)
% end
% 
% plot(ExpectedValues1./TotalNuclei,ExpectedValues2./TotalNuclei,'r-','LineWidth',2)
% %hold off
% 
% xlabel('nuclei with one spot')
% ylabel('nuclei with two spots')
% % ylim([0 1])
% % xlim([0 1])
% set(gca,'FontSize',18)

end

