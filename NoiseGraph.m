
dA = 0.7
dB = 0.1

close all
hold on
plot([0 2],[0 2],'k-')
plot([0 2],[1+dA 1+dA],'r--')
plot([0 2],[1+dB 1+dB],'k--')
plot([1+dB 1+dB],[0 2],'k--')
plot([1+dA 1+dA],[0 2],'r--')
plot([1+dB 1+(dA+dB)/2],[1+dA 1+(dA+dB)/2],'g-','LineWidth',2)
plot([1 1+dB],[1 1+dA],'b-','LineWidth',2)
plot([1 1+(dA+dB)/2],[1 1+(dA+dB)/2],'k-','LineWidth',2)
plot(1+dB,1+dA,'ro','MarkerFaceColor','r','MarkerSize',7)
plot(1,1,'ko','MarkerFaceColor','k','MarkerSize',7)
plot(1+(dA+dB)/2,1+(dA+dB)/2,'bo','MarkerFaceColor','b','MarkerSize',7)
hold off
ylabel('A (normalized to the mean)')
xlabel('B (normalized to the mean)')
legend('A=B','1+dA','1+dB','','','d_{int}','d_{tot}','d_{corr}',...
    'single cell','population mean','projection onto A=B')