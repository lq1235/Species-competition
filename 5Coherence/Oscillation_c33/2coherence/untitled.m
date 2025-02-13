clc
clear
a=xlsread("data.xlsx");
plot(a(:,1),a(:,2),"b-o",'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 10)
xlabel('\fontsize{27} C_{33}');
ylabel('\fontsize{27} Coherence');


ax = gca();
% ax.YRuler.Exponent = -1;
ax.XRuler.Exponent = -1;
set(gca,'xtick',0.125:0.01:0.165)
xlim([0.122441113044196,0.142591113044196])
ylim([0.622321631524838,0.864321631524837])
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,0 0.00025])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'xtick',0.125:0.005:0.165)
