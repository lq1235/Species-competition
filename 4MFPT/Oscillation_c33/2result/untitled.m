clc
clear
a=load('MFPT.txt');
%MFPT
h=figure(1)
plot(a(:,1), a(:,2), '^--', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,-0.0035 0.04])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.YRuler.Exponent = 4;
ax.XRuler.Exponent = -2;
% set(gca,'xtick',0.056:0.0005:0.059)
xlabel("\fontsize{27} C_{33}");
ylabel("\fontsize{27} MFPT*");
xlim([0.122058117028903,0.144019617028903])
ylim([-179481.3589780138,991798.6410219865])

ax = gca();
ax.YRuler.Exponent = 5;
ax.XRuler.Exponent = -1;

h=figure(2)
plot(a(:,1), log(a(:,2)), '^--', 'Color', 'g', 'LineWidth', 1, 'MarkerFaceColor', 'g', 'MarkerSize', 10)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,-0.0035 0.04])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
% set(gca,'xtick',0.056:0.0005:0.059)
xlabel("\fontsize{27} C_{33}");
ylabel("\fontsize{27} MFPT*");
xlim([0.12142700129493,0.14338850129493])
ylim([1.873837277752573,17.117422203519762])

ax = gca();
ax.YRuler.Exponent = 0;
ax.XRuler.Exponent = -1;

