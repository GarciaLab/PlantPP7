function CalculateNoiseComponents(Array)

% calculates the total, intrinsic and extrinsic noise based on the
% treatment by Elowitz et al. 2002

% ARRAY is a a x n matrix where a = number of alleles per cell (hard coded
% to be 2 for now) and n is the number of cells. It contains the expression level (rows) of
% each cell (columns).

FirstAlleleData = Array(1,:);
SecondAlleleData = Array(2,:);

% NoiseTotal = 1/2 * (mean(FirstAlleleData.^2) + mean(SecondAlleleData.^2));
% NoiseIntrinsic = 1/2 * mean((FirstAlleleData - SecondAlleleData).^2);
% NoiseCorrelated = mean(FirstAlleleData.*SecondAlleleData);

% title({Prefix,['\eta_{tot}^2 =' num2str(NoiseTotal) ' \eta_{int}^2 =' num2str(NoiseIntrinsic) ...
%     ' \eta_{corr}^2 =' num2str(NoiseCorrelated)]})
figure
plot(FirstAlleleData,SecondAlleleData,'o','Color','k','MarkerFaceColor','r','MarkerSize',13)
hold on
plot([0 max(FirstAlleleData,SecondAlleleData)*1.1],[0 max(FirstAlleleData,SecondAlleleData)*1.1],'k-')
xlabel('integrated fluorescence allele 1 (normalized to mean of all alleles)')
ylabel('integrated fluorescence allele 2 (normalized to mean of all alleles)')

% calculate deltaA and deltaB, the distance of each normalized point to the
% mean of 1.
FirstAlleleData = 1 - FirstAlleleData;
SecondAlleleData = 1- SecondAlleleData;
%% bootstrap the error in noise components to get errobars
NoiseTotal = 1/2 * (mean(FirstAlleleData.^2) + mean(SecondAlleleData.^2));
NoiseIntrinsic = 1/2 * mean((FirstAlleleData - SecondAlleleData).^2);
NoiseCorrelated = mean(FirstAlleleData.*SecondAlleleData);

nSamples = 10000; %how many times we're bootstrapping
%these are the function we are bootstrapping, the formulas for each noise component
TotObsP = @(x,y) 1/2 * (mean(x.^2) + mean(y.^2)); 
IntObsP = @(x,y) 1/2 * mean((x-y).^2);
CorrObsP = @(x,y) mean(x.*y);

% bootstrapping the observed total noise
[bootTotObs,~] = bootstrp(nSamples,TotObsP,FirstAlleleData,SecondAlleleData); 
bootstrappedMeanTotalNoise = mean(bootTotObs);
bootstrappedStdTotalNoise = std(bootTotObs);

% bootstrapping the observed intrinsic noise
[bootIntObs,~] = bootstrp(nSamples,IntObsP,FirstAlleleData,SecondAlleleData); 
bootstrappedMeanIntrinsicNoise = mean(bootIntObs);
bootstrappedStdIntrinsicNoise = std(bootIntObs);

% bootstrapping the observed correlated noise
[bootCorrObs,~] = bootstrp(nSamples,CorrObsP,FirstAlleleData,SecondAlleleData); 
bootstrappedMeanCorrNoise = mean(bootCorrObs);
bootstrappedStdCorrNoise = std(bootCorrObs);

figure
hold on
errorbar(1,bootstrappedMeanTotalNoise,bootstrappedStdTotalNoise,'bo','LineWidth',3,...
    'CapSize',0)
errorbar(2,bootstrappedMeanIntrinsicNoise,bootstrappedStdIntrinsicNoise,'ro','LineWidth',3,...
    'CapSize',0)
errorbar(3,bootstrappedMeanCorrNoise,bootstrappedStdCorrNoise,'go','LineWidth',3,...
    'CapSize',0)
hold off
ylim([0 1])
xlim([0.8 3.2])
legend('\eta_{tot}^2','\eta_{int}^2','\eta_{ext}^2')


