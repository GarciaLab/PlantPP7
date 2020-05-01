function plotSanityCheck(Struct,Field1,Field2,Field3)

%this function plots a scatter of Field1 x Field2 versus Field3
% use as plotSanityCheck(alignedDatasetsStruct,'InstFractionON','MeanFluoOn','MeanFluoAll')


Field1All = plotAllPrefixes(Struct,Field1); 
Field2All = plotAllPrefixes(Struct,Field2); 
Field3All = plotAllPrefixes(Struct,Field3);
%close all %close the figures that the plotAllPrefixes function generates

figure
Expectation = Field1All.*Field2All;
hold on
scatter(Field3All(:),Expectation(:),200,'b','filled','MarkerFaceAlpha',0.15);
plot(Field3All(:),Field3All(:),'k-','LineWidth',1.5)
hold off
xlabel(Field3)
ylabel([Field1 ' x ' Field3])
legend('data','y=x')
title('Sanity Check')