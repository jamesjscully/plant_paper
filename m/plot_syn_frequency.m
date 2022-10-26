figure(1)
clf 
size(fr1)
size(avs)
plot (fr1,1.6*avs(1:end),'.-','MarkerSize',10,'Color',[0.5 0.5 0.5])
hold on
plot (fr1,5*avs1(1:end),'.-','MarkerSize',10,'Color',[0.1 0.7 0.7])
hold on
plot (fr1,1.2*avs2(1:end),'.-','MarkerSize',10,'Color',[0.9 0. 0.0])
hold on
plot (fr1,1.6*avs3(1:end),'.-','MarkerSize',10,'Color',[0 0. 0.0])
hold on
 ylim([0 1]) 
 xlim([0 11]) 
   ax = gca;
 ax.FontSize = 14;
 xlabel('frequency [Hz]','Fontsize',18),ylabel('neurotrans-rate','Fontsize',18)
