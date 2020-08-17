function [X,Y,TwoSpots2SDs] = MultinomialTestOfTwoSpots(FractionOneSpot,FractionTwoSpots,TotalNuclei,bootstrap)

%% Binomial distribution - analytical approach

C = [0 0 1]; %color for probability shades
w = [1 1 1]; %to dilute color C
n = TotalNuclei; %number of independent trials

%p = 0.5;
valuesForP = [0:0.005:1];
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
    
    TwoSpots2SDs(counter) = sqrt(Var2)*2;
    Yup(counter) = E2+sqrt(Var2)*2;
    Ydown(counter) = E2-sqrt(Var2)*2;
%     X3(counter) = E1+sqrt(Var2)*2;
%     X4(counter) = E1-sqrt(Var2)*2;
    ExpectedValues1(counter) = E1;
    ExpectedValues2(counter) = E2;
    X(counter) = E1;
    Y(counter) = E2;
    
    counter = counter+1;

end

if bootstrap
    % plot data bootstrapping it
    [boostrpMean1Spot,boostrpErr1Spot] = bootstrapNSpots(n,ceil(FractionOneSpot*n));
    [boostrpMean2Spot,boostrpErr2Spot] = bootstrapNSpots(n,ceil(FractionTwoSpots*n));

    errorbar(boostrpMean1Spot,boostrpMean2Spot,boostrpErr1Spot,'horizontal','k',...
        'CapSize',0,'LineWidth',2)
    errorbar(boostrpMean1Spot,boostrpMean2Spot,boostrpErr2Spot,'k',...
        'CapSize',0,'LineWidth',2)
    plot(boostrpMean1Spot,boostrpMean2Spot,'ko','Color','k','MarkerSize',7,...
        'MarkerFaceColor','w','LineWidth',2)
    
else
    plot(FractionOneSpot*n,FractionTwoSpots*n,'ko','Color','k','MarkerFaceColor','k',...
        'MarkerSize',5,'LineWidth',2)
end

plot(ExpectedValues1,ExpectedValues2,'r-','LineWidth',2)
plot(X,Yup,'-','Color','b')
plot(X,Ydown,'-','Color','b')
xlabel('nuclei with one spot')
ylabel('nuclei with two spots')
set(gca,'FontSize',18)

end

