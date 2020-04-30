function AUperGFP = AbsoluteCalibration(FluorescenceVector60,FluorescenceVector120,...
    FluoFactor60,FluoFactor120)


figure(2)
Xvalues = [60 120]; 

[Fitted60mersMean Fitted60mersSD]  = fitsToNanocages(FluorescenceVector60,FluoFactor60); %fits data to a gaussian and returns mean and sd
[Fitted120mersMean Fitted120mersSD]= fitsToNanocages(FluorescenceVector120,FluoFactor120);

Yvalues = [Fitted60mersMean Fitted120mersMean]; %these are the means corresponding to each X value
SDs = [Fitted60mersSD Fitted120mersSD]; %these are the standard deviations of each mean

%now standard errors for plotting only
SEMs = [std(FluorescenceVector60./FluoFactor60)/sqrt(length(FluorescenceVector60)) ...
    std(FluorescenceVector120./FluoFactor120)/sqrt(length(FluorescenceVector120))];




%%

% hold on
% SlopeThroughOrigin = [0 X]'\[0 Y]';
% FitY = SlopeThroughOrigin .* [0 X 130];
% plot([0 X 130],FitY,'r-','LineWidth',2)
% errorbar(X(1),Y(1),SEM(1),'bo','LineWidth',2,'CapSize',0,'MarkerFaceColor','b')
% errorbar(X(2),Y(2),SEM(2),'ko','LineWidth',2,'CapSize',0,'MarkerFaceColor','k')
% xlim([0 130])
% xticks([60 120])
% legend(['y=' num2str(round(SlopeThroughOrigin,3)) 'x'],'60 GFP \pm sd','120 GFP  \pm sd')
% hold off
% ylabel('fluorescence (AU)')
% xlabel('number of GFP molecules')
% set(gca,'FontSize', 18,'FontWeight','Bold')
% title({'Absolute GFP calibration \sim ',[num2str(round(1/round(SlopeThroughOrigin,3),2)) ' GFPs / AU']})
% 
% %% plot linear fit without forcing through origin
% figure(3)
% 
% Data120mer = AllFixedAreaIntensityVector_120./3;
% Data60mer = AllFixedAreaIntensityVector_60./5;
% 
% X = [0 60*ones(size(Data60mer)) 120*ones(size(Data120mer))];
% Y = [0 Data60mer Data120mer];
% [p,S] = polyfit(X,Y,1); 
% [y_fit,delta] = polyval(p,X,S);
% 
% plot(X,y_fit,'r-')
% hold on
% plot(X,y_fit+2*delta,'m--',X,y_fit-2*delta,'m--')
% 
% 
% 
% %% regression
% %problem: there's no data for 0 GFPs so the regression doesn't care about
% %the 0. I can fix this by creating 'fake data' for 0GFPs, with a fluorescence
% %of 0 a.u. Because by definition 0GFPs = 0 a.u we'll create a ridiculous
% %number of '0GFP nanocage' measurements.
% 
% 
% x2 = [0*ones(1,5000) 60*ones(size(Data60mer)) 120*ones(size(Data120mer))]';
% Y2 = [0*ones(1,5000) Data60mer Data120mer]';
% 
% X2 = [ones(size(x2)) x2 ];% the padding with ones is needed for the thing to run
% [b,bint,r,rint,stats] = regress(Y2,X2) 
% slope = b(2);
% minslope = bint(2,2) ;
% maxslope = bint (2,1);
% Yintersect = b(1);
% 
% Xs = [60 120];
% Ys = [Guess60(2) Guess120(2)];
% Es = [Guess60(3) Guess120(3)];
% 
% hold on
% plot(x2,Yintersect+x2*slope,'r-','LineWidth',3)
% plot(x2,Yintersect+x2*minslope,'r:')
% plot(x2,Yintersect+x2*maxslope,'r:')
% errorbar(Xs(1),Ys(1),Es(1),'bo','LineWidth',2,'CapSize',0,'MarkerFaceColor','b')
% errorbar(Xs(2),Ys(2),Es(2),'ko','LineWidth',2,'CapSize',0,'MarkerFaceColor','k')
% xlim([0 130])
% 
% 
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


%%

% X = [zeros(1,50000000) ones(1,length(FluorescenceVector60))*60 ones(1,length(FluorescenceVector120))*120]; %these are the X values
% Y = [zeros(1,50000000) FluorescenceVector60/FluoFactor60 FluorescenceVector120/FluoFactor120]; %and the Y corresponding to each X value
% % note that we're assuming that zero GFPs have a fluorescence of zero a.u
% 
% slope = X(:)\Y(:); % linear regression throug origin
% AUperGFP = slope;
% yFit = X(:)*slope; % calculate fitted line
% 
% % calculate R2          
% f = yFit;
% Bbar = mean(X');
% SStot = sum((X' - Bbar).^2);
% SSreg = sum((f - Bbar).^2);
% SSres = sum((X' - f).^2);
% R2 = 1 - SSres/SStot;
% R = corrcoef(Y,X');
% Rsq = R(1,2).^2
% 
% %% plot results
% 
% figure
% hold on
% % the data
% errorbar(Xvalues,Yvalues,SEMs,'o','LineWidth',2,'CapSize',0)
% % the fit
% plot([0 Xvalues],[0 Xvalues].*slope,'k-')
% 
% hold off


%% regression

%problem: there's no data for 0 GFPs so the regression doesn't care about
%the 0. I can fix this by creating 'fake data' for 0GFPs, with a fluorescence
%of 0 a.u. Because by definition 0GFPs = 0 a.u we'll create a ridiculous
%number of '0GFP nanocage' measurements.


x2 = [0*ones(1,500000) ones(size(FluorescenceVector60))*60 ones(size(FluorescenceVector120))*120]'; %these are the X values
Y2 = [0*ones(1,500000) FluorescenceVector60/FluoFactor60 FluorescenceVector120/FluoFactor120]'; %and the Y corresponding to each X value

X2 = [ones(size(x2)) x2 ];% the padding with ones is needed for the thing to run
[b,bint,r,rint,stats] = regress(Y2,X2);
Rsqr = round(stats(1),2);
slope = round(b(2),3);
AUperGFP = slope;
minslope = bint(2,2) ;
maxslope = bint (2,1);
Yintercept = b(1);

Xs = Xvalues;
Ys = Yvalues;
Es = SDs./sqrt([length(FluorescenceVector60) length(FluorescenceVector120)]);

figure
hold on
plot(x2,Yintercept+x2*slope,'r-','LineWidth',3)
errorbar(Xs(1),Ys(1),Es(1),'bo','LineWidth',2,'CapSize',0,'MarkerFaceColor','b')
errorbar(Xs(2),Ys(2),Es(2),'ko','LineWidth',2,'CapSize',0,'MarkerFaceColor','k')
xlim([0 130])
ylabel('fluorescence (a.u)')
xlabel('number of GFP molecules')
legend(['linear regression Y=' num2str(Yintercept) '+' num2str(slope) 'X ; R^2=' num2str(Rsqr)],...
    '60GFP +- SEM','120GFP +- SEM')



end
