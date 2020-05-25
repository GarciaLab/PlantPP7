

%% predict F1 and F2 from F0
F0 = linspace(0,1,100);
p0 = 1-sqrt(F0);
F1 = (2.*p0).*(1-p0);
F2 = p0.^2;

figure(1)
hold on
yyaxis left
plot(F0,F1)
ylabel('F1')
ylim([0 1])
yyaxis right
plot(F0,F2)
ylabel('F2')
hold off
xlabel('F0')


%% predict F0 and F2 from F1

F1 = linspace(0,.5,100);
p1 = 1/2.*(1-sqrt(1-2.*F1));
F0 = (1-p1).^2;
F2 = p1.^2;

figure(2)
hold on
yyaxis left
plot(F1,F0)
ylabel('F0')
ylim([0 1])
yyaxis right
plot(F1,F2)
ylabel('F2')
hold off
xlabel('F1')

%% predict F0 and F1 from F2

F2 = linspace(0,1,100);
F0 = (1-sqrt(F2)).^2;
F1 = 2.*sqrt(F2) - 2.*F2;

figure(3)
hold on
yyaxis left
plot(F2,F0)
ylabel('F0')
ylim([0 1])
yyaxis right
plot(F2,F1)
ylabel('F1')
ylim([0 1])
hold off
xlabel('F2')


%%
theArbitraryNumber = 2;
valuesForP = 0:10^-theArbitraryNumber:1;
% find where 0.5 occurs because F1 can't go further than that
[~,idx] = find(valuesForP==0.5);
% create a matrix to store results
% the first dimension will be F0, the second F1 and the entries will
% reflect the value of F2
FractionsMatrix = zeros(length(valuesForP),idx);
valuesForF0 = (1-valuesForP).^2;

valuesForF1 = [];
valuesForF2 = [];
for F0 = valuesForF0
    
    F1 = 2*sqrt(F0) - 2*F0; %prediction of F1 from F0
    [~,F1idx] = min(abs(valuesForP-F1));
    F2 = F0 - 2*sqrt(F0) + 1; %prediction of F2 from F0
    %[~,F2idx] = min(abs(valuesForP-F2));
    [~,F0idx] = find(valuesForF0==F0);
    FractionsMatrix(F0idx,F1idx) = F2;
    valuesForF1 = [valuesForF1 F1];
    valuesForF2 = [valuesForF2 F2];
end


%%
figure

% p3data = plot3(0.5,0.3,0.2,'bo');
% p3data.MarkerFaceColor = 'b';

p3 = plot3(valuesForF0,valuesForF1,valuesForF2,'r-o');
p3.LineWidth = 2;
hold on
p3data = plot3(0.5,0.3,0.2,'bo');
p3data.MarkerFaceColor = 'b';
hold off
xlabel('F_0')
ylabel('F_1')
zlabel('F_2')
grid on

%% Probability distribution 
% of obtaining two heads Nhh times and one head
%Nht times in Ntot flips

Ntot = 559;
NhhValues = 0:1:Ntot;
hold on
p = 0.5; %this is the probability of getting heads in one coin flip
phh = p^2;
pht = 2*p*(1-p);
counterNhh = 1;
Matrix = zeros(length(NhhValues),length(NhhValues));
for Nhh = NhhValues  
    P_HH = phh^(Nhh) * (1-phh)^(Ntot-Nhh) * nchoosek(Ntot,Nhh); 
    counterNht = 1;
    for Nht = 0:Ntot-Nhh
        P_HT = pht^(Nht) * (1-pht)^(Ntot-Nht) * nchoosek(Ntot,Nht);
        P_both = P_HH * P_HT;
        Matrix(counterNhh,counterNht) = P_both;       
        counterNht = counterNht+1;
    end
    counterNhh = counterNhh+1;
end

Matrix = Matrix./max(Matrix(:));
contourf(Matrix,'ShowText','on')
ylabel('two spots')
xlabel('one spot')
zlabel('probability')
% xlim([40 70])
% ylim([10 50 ])
%% Probability distribution 
% of obtaining two heads Nhh times and one head
%Nht times in Ntot flips...for multiple values of p

Ntot = 559; %upper bound: all nuclei can potentially have a spot
%Ntot = ceil(0.1984*Ntot + 0.0754*Ntot); %lower bound: all the nuclei that could have a spot have at least one spot
NhhValues = 1:Ntot;
BigMatrix = zeros(length(NhhValues),length(NhhValues));
BigCountMatrix = BigMatrix;

hold on
for p = 0.1:0.5:0.9 %this is the probability of getting heads in one coin flip
    phh = p^2;
    pht = 2*p*(1-p);
    counterNhh = 1;
    Matrix = zeros(length(NhhValues),length(NhhValues));
    for Nhh = NhhValues 
        P_HH = phh^(Nhh) * (1-phh)^(Ntot-Nhh) * nchoosek(Ntot,Nhh); 
        counterNht = 1;
        for Nht = 1:Ntot-Nhh
            P_HT = pht^(Nht) * (1-pht)^(Ntot-Nht) * nchoosek(Ntot,Nht);
            P_both = P_HH * P_HT;
            Matrix(counterNhh,counterNht) = P_both;       
            counterNht = counterNht+1;
        end
        counterNhh = counterNhh+1;
    end
    
    NormMatrix = Matrix./max(Matrix(:)); %normalize by the mean
    BinNormMatrix = NormMatrix > 0.001; 
    BigCountMatrix = BigCountMatrix + BinNormMatrix;
    BigMatrix = BigMatrix + NormMatrix;
end

BigCountMatrix(BigCountMatrix==0)=nan;
close all
figure
hold on
contourf(BigMatrix./BigCountMatrix,'ShowText','on')
ylabel('two spots (% nuclei)')
xlabel('one spot (% nuclei)')
plot(0.1984*Ntot,0.0754*Ntot,'ko','MarkerSize',10,'MarkerFaceColor','w')
% plot(0.18*100,0.05*100,'bo','MarkerSize',10,'MarkerFaceColor','w')
% plot(0.1*100,0.07*100,'bo','MarkerSize',10,'MarkerFaceColor','w')
% plot(0.18*100,0.18*100,'ro','MarkerSize',10,'MarkerFaceColor','w')

colormap viridis
set(gca,'FontSize',18)
hold off

%% Binomial distribution - analytical approach



n = 400; %number of independent trials
%p = 0.5;
valuesForP = [0:0.1:1];
Expected1spot = nan(1,length(valuesForP));
Expected2spot = nan(1,length(valuesForP));
Error1spot = nan(1,length(valuesForP));
Error2spot = nan(1,length(valuesForP));

counter = 1;
for p = valuesForP
    % x0 = ; %random variable: number of cells with zero spots
    % x1 = ; %random variable: number of cells with 1 spot
    % x2 = ; %random variable: number of cells with 2 spots
    p0 = (1-p)^2; %probability of x0
    p1 = 2*p*(1-p); %probability of x1
    p2 = p^2; %probability of x2

    %E0 = n * p0; %expected value of cells with zero spots
    E1 = n * p1; %expected value of cells with one spots
    E2 = n * p2; %expected value of cells with two spots
    Expected1spot(counter) = E1;
    Expected2spot(counter) = E2;

    %Var0 = n * p0 * (1-p0); %variance in the number of cells with zero spots
    Var1 = n * p1 * (1-p1); %variance in the number of cells with one spots
    Var2 = n * p2 * (1-p2); %variance in the number of cells with two spots
    Error1spot(counter) = sqrt(Var1);
    Error2spot(counter) = sqrt(Var2);
    %CoVar1_2 = -n*p1*p2; % covariance in the number of cells with one spot and the number of cells with two spots   
    counter = counter+1;

end

hold on
errorbar(Expected1spot,Expected2spot,2.*Error1spot,'CapSize',0,'LineStyle','none',...
    'LineWidth',2)
errorbar(Expected1spot,Expected2spot,2.*Error2spot,'horizontal','CapSize',0,...
    'LineStyle','none','LineWidth',2)
plot(0.1984*Ntot,0.0754*Ntot,'ko','MarkerSize',10,'MarkerFaceColor','w')
hold off
ylim([0 Ntot])
xlim([0 Ntot])
xlabel('nuclei with one spot')
ylabel('nuclei with two spots')