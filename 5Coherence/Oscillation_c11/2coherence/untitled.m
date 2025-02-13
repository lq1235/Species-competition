clc
clear
a=xlsread("data.xlsx");
plot(a(:,1),a(:,2),"b-o",'LineWidth', 1, 'MarkerFaceColor', 'blue', 'MarkerSize', 10)

xlabel('\fontsize{27} C_{11}');
ylabel('\fontsize{27} Coherence');


ax = gca();
% ax.YRuler.Exponent = -1;
ax.XRuler.Exponent = -2;
% set(gca,'xtick',0.125:0.01:0.165)
xlim([0.055818283965362,0.058])
ylim([0.892718960484473,1.012508960484473])
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,0 0.00025])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'xtick',0.0550:0.0005:0.0615)