function OutputData = plotRTqPCRResults(path);

load(path)

%Combine all three replicates in one figure

YHSP101_means_all = mean([YHSP101_means_3,YHSP101_means_4,YHSP101_means_5],2);
PP7LUC_means_all = mean([PP7LUC_means_3,PP7LUC_means_4,PP7LUC_means_5],2);



%propagate error!
YHSP101_means = [YHSP101_means_3,YHSP101_means_4,YHSP101_means_5];
YHSP101_SE_all = std(YHSP101_means')./sqrt(3);

PP7LUC_means = [PP7LUC_means_3,PP7LUC_means_4,PP7LUC_means_5];
PP7LUC_SE_all = std(PP7LUC_means')./sqrt(3);

figure
hold on
errorbar(time,YHSP101_means_all,YHSP101_SE_all,'b','LineWidth',2,'CapSize',0)
errorbar(time,PP7LUC_means_all,PP7LUC_SE_all,'r','LineWidth',2,'CapSize',0)
hold off
xlabel('time at 38C (min)')
ylabel('Expression level (Act2 units)')
title('Mean of replicates 3,4 and 5')
%set(gca,'fontsize', 18,'FontWeight','Bold')
legend ('mean HSP101 mRNA and SEM','mean reporter mRNA and SEM')
xlim([0 65])


%% regression
LabelsForPoints = {'0','5','10','15','30','60'};
transpYHSP101_means = YHSP101_means'; %transpose before converting into a vector
Xdata = transpYHSP101_means(:); 
transpPP7_means = PP7LUC_means'; %transpose before converting into a vector
Ydata = transpPP7_means(:); 

regX = [ones(size(Xdata)) Xdata];% the padding with ones is needed for the thing to run
[b,bint,r,rint,stats] = regress(Ydata,regX);
slope = b(2);
minslope = bint(2,2) ; % top of 95% confidence interval
maxslope = bint (2,1);% bottom of 95% confidence interval
Yintersect = b(1);

figure
hold on
% plot the fit results
plot(Xdata,Yintersect+Xdata*slope,'k-','LineWidth',1)
%plot(Xdata,Yintersect+Xdata*minslope,'-','Color',[.5 .5 .5])
%plot(Xdata,Yintersect+Xdata*maxslope,'-','Color',[.5 .5 .5])
%now the data

% double error bar, endogenous HSP101 on the x, PP7 on the y
errorbar(YHSP101_means_all,PP7LUC_means_all,YHSP101_SE_all,'b.','horizontal','CapSize',0,'LineWidth',2)
errorbar(YHSP101_means_all,PP7LUC_means_all,PP7LUC_SE_all,'r.','CapSize',0,'LineWidth',2)
% plot the center of the double errorbars
plot(YHSP101_means_all,PP7LUC_means_all,'ko','MarkerSize',8,'MarkerFaceColor','w')
%label the points with the time under induction
labelpoints(YHSP101_means_all,PP7LUC_means_all,LabelsForPoints,'NW',1,0.6)

hold off
xlabel('endogenous HSP101 mRNA abundance (Act2 units)')
ylabel('transgene mRNA abundance (Act2 units)')
legend(['linear regression, slope=' num2str(slope) ': R^2=' num2str(stats(1))],...
    'HSP101 mean +- SEM','transgene +- SEM')
title('mRNA abundance, transgene vs endogenous HSP101')
Limits = xlim;
ylim(Limits)

%% now normalizing to t=60 min
%divide CTs by the mean value across replicates at 60min
norm_YHSP101_means = YHSP101_means./mean(YHSP101_means(end,:));
YHSP101_means_All_Norm = mean(norm_YHSP101_means'); 
YHSP101_SE_All_Norm = std(norm_YHSP101_means')/sqrt(size(norm_YHSP101_means,1));
norm_PP7LUC_means = PP7LUC_means./mean(PP7LUC_means(end,:));
PP7_means_All_Norm = mean(norm_PP7LUC_means');
PP7_SE_All_Norm = std(norm_PP7LUC_means')/sqrt(size(norm_PP7LUC_means,1));

transpYHSP101_means = norm_YHSP101_means'; %transpose before converting into a vector
Xdata = transpYHSP101_means(:); 
transpPP7_means = norm_PP7LUC_means'; %transpose before converting into a vector
Ydata = transpPP7_means(:); 

regX = [ones(size(Xdata)) Xdata];% the padding with ones is needed for the thing to run
[b,bint,r,rint,stats] = regress(Ydata,regX);
slope = b(2);
minslope = bint(2,2) ; % top of 95% confidence interval
maxslope = bint (2,1);% bottom of 95% confidence interval
Yintersect = b(1);

figure
hold on
% plot the fit results
plot(Xdata,Yintersect+Xdata*slope,'k-','LineWidth',1)

% double error bar, endogenous HSP101 on the x, PP7 on the y
errorbar(YHSP101_means_All_Norm,PP7_means_All_Norm,YHSP101_SE_All_Norm,'b.','horizontal','CapSize',0,'LineWidth',2)
errorbar(YHSP101_means_All_Norm,PP7_means_All_Norm,PP7_SE_All_Norm,'r.','CapSize',0,'LineWidth',2)
% plot the center of the double errorbars
plot(YHSP101_means_All_Norm,PP7_means_All_Norm,'ko','MarkerSize',8,'MarkerFaceColor','w')
labelpoints(YHSP101_means_All_Norm,PP7_means_All_Norm,LabelsForPoints,'NW',1,0.6)

hold off
xlabel('endogenous HSP101 mRNA abundance (normalized to 60 min)')
ylabel('transgene mRNA abundance (normalized to 60 min)')
legend(['linear regression, slope=' num2str(slope) ': R^2=' num2str(stats(1))],...
    'HSP101 mean +- SEM','transgene +- SEM','Location','NorthWest')
title('mRNA abundance, transgene vs endogenous HSP101 - normalized to t=60min')
Limits = xlim;
ylim(Limits)

%% fit reporter mRNA abundance (mean across reps) to mRNA prodution-synthesis model
Xdata = [0 5 10 15 30 60];
Ydata = PP7LUC_means_all;
% the function: the solution to dM/dt = r - gamma*M which is
% M(t) = r/gamma + c * e^(-C*t) where C is a constant determined by the
% initial number of mRNAs. if the initial number is 0 then c = -r/gamma

mRNAfun = @(z) (((-z(1)/z(2)) * exp(-z(1)*Xdata)) + (z(1)/z(2))) - Ydata;
% z(1) is the transcription rate r, z(2) is the degradation rate gamma

z0 = [0.01,0.000000001]; %initial guesses
lb = [0.02,0.000001]; % lower bounds
ub = [0.1,0.006]; % define upper bounds for each parameter, note that I'm forcing the baseline to be 0
Guess = lsqnonlin(mRNAfun,z0,lb,ub); %fit by minizing the subtraction between actual Y data and the fit
moreContinuousXData = linspace(min(Xdata),max(Xdata)*1.2,100);
FittedResult = (-Guess(1)/Guess(2))*exp(-Guess(2)*moreContinuousXData) + (Guess(1)/Guess(2));
plot(moreContinuousXData,FittedResult,'b','LineWidth',2);
hold on
%plot(Xdata,PP7LUC_means_all,'ro')
errorbar(time,PP7LUC_means_all,PP7LUC_SE_all,'ro','LineWidth',2,'CapSize',0)
mRNALifeTime = log(2)/Guess(2); %in minutes
legend(['r =' num2str(Guess(1)) '; \gamma =' num2str(Guess(2))])
title(['mRNA lifetime = ' num2str(mRNALifeTime) 'min'])
%hold off

%% function output
OutputData = norm_PP7LUC_means';
end