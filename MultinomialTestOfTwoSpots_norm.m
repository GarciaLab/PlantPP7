function MultinomialTestOfTwoSpots_norm(FractionOneSpot,FractionTwoSpots,TotalNuclei,bootstrap)

%% Binomial distribution - analytical approach

C = [0 0 1]; %color for probability shades
w = [1 1 1]; %to dilute color C
n = TotalNuclei; %number of independent trials

%p = 0.5;
valuesForP = [0.01:0.005:0.99];
Expected1spot = nan(1,length(valuesForP));
Expected2spot = nan(1,length(valuesForP));
Error1spot = nan(1,length(valuesForP));
Error2spot = nan(1,length(valuesForP));

counter = 1;
Yup= [];
Ydown= [];
X =[];
Y2 =[];
for p = valuesForP
    p1 = 2*p*(1-p); %probability of x1
    p2 = p^2; %probability of x2
    E1 = n * p1; %expected value of cells with one spots
    E2 = n * p2; %expected value of cells with two spots
    Expected1spot(counter) = E1;
    Expected2spot(counter) = E2;
    Var1 = n * p1 * (1-p1); %variance in the number of cells with one spots
    Var2 = n * p2 * (1-p2); %variance in the number of cells with two spots
    Error1spot(counter) = sqrt(Var1);
    Error2spot(counter) = sqrt(Var2);

    Yup(counter) = E2+sqrt(Var2)*2;
    Ydown(counter) = E2-sqrt(Var2)*2;
%     X3(counter) = E1+sqrt(Var2)*2;
%     X4(counter) = E1-sqrt(Var2)*2;
    ExpectedValues1(counter) = E1;
    ExpectedValues2(counter) = E2;
    X(counter) = E1;
    
    counter = counter+1;

end

if bootstrap
    % plot data bootstrapping it
    [boostrpMeanFrac1Spot,boostrpErrFrac1Spot] = bootstrapFracNSpots(n,ceil(FractionOneSpot*n));
    [boostrpMeanFrac2Spot,boostrpErrFrac2Spot] = bootstrapFracNSpots(n,ceil(FractionTwoSpots*n));

    errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac1Spot,'horizontal','k',...
        'CapSize',0,'LineWidth',2)
    errorbar(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,boostrpErrFrac2Spot,'k',...
        'CapSize',0,'LineWidth',2)
    plot(boostrpMeanFrac1Spot,boostrpMeanFrac2Spot,'ko','Color','k','MarkerSize',7,...
        'MarkerFaceColor','w','LineWidth',2)
    
else
    plot(FractionOneSpot,FractionTwoSpots,'ko','Color','k','MarkerFaceColor','k',...
        'MarkerSize',5,'LineWidth',2)
end

plot(ExpectedValues1./TotalNuclei,ExpectedValues2./TotalNuclei,'r-','LineWidth',2)
plot(X./TotalNuclei,Yup./TotalNuclei,'-','Color','b')
plot(X./TotalNuclei,Ydown./TotalNuclei,'-','Color','b')
%hold off

xlabel('nuclei with one spot')
ylabel('nuclei with two spots')
% ylim([0 1])
% xlim([0 1])
set(gca,'FontSize',18)

end

