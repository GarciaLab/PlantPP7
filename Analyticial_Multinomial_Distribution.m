

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

F1 = figure;
hold on

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
    counter = counter+1;
    
    % plot a shaded area of 1,2 and 3 standard deviations
    drawellipse('Center',[E1,E2],'SemiAxes',4*[sqrt(Var1),sqrt(Var2)],'Color',(C+16*w)/17,...
        'handlevisibility','off','InteractionsAllowed','none','LineWidth',0.00001,...
        'FaceAlpha',1);
end

MinusTwoSDLine = [];
PlusTwoSDLine = [];
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
    counter = counter+1;
    
    % plot a shaded area of 1,2 and 3 standard deviations
    drawellipse('Center',[E1,E2],'SemiAxes',3*[sqrt(Var1),sqrt(Var2)],'Color',(C+6*w)/7,...
        'handlevisibility','off','InteractionsAllowed','none','LineWidth',0.00001,...
        'FaceAlpha',1);
end


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
    counter = counter+1;
    
    % plot a shaded area of 1,2 and 3 standard deviations
    drawellipse('Center',[E1,E2],'SemiAxes',2*[sqrt(Var1),sqrt(Var2)],'Color',(C+3*w)/4,...
        'handlevisibility','off','InteractionsAllowed','none','LineWidth',0.00001,...
        'FaceAlpha',1);
end

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
    
    % plot a shaded area of 1,2 and 3 standard deviations
    drawellipse('Center',[E1,E2],'SemiAxes',[sqrt(Var1),sqrt(Var2)],'Color',(C+2*w)/3,...
        'handlevisibility','off','InteractionsAllowed','none','LineWidth',0.00001,...
        'FaceAlpha',1);

end


hold off

ylim([0 n])
xlim([0 n])
xlabel('nuclei with one spot')
ylabel('nuclei with two spots')
set(gca,'FontSize',18)

end

