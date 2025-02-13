clc
clear
a=load('MFPT.txt');
%MFPT
% h=figure(1)
% plot(a(:,1),a(:,2),'ko-','LineWidth',1,'markersize',10)
% hold on
% plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
% hold on

% axis([-0.5 3.5 ,-1*10^(-3) 13*10^(-3)])
% set(gca,'LineWidth',1.2,'Fontsize',27)
% set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% xlabel("\fontsize{27}C_{15}");
% ylabel("\fontsize{27}T_{MFPT}");
% plot([0.7 0.7],[-0.5*10^(5) 9*10^(5)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-0.5*10^(5) 9*10^(5)],'b--','LineWidth',1)
% axis([-0.2 3.2,-0.5*10^(5) 9*10^(5)])
% xlim([0.001 0.0025])


%1/MFPT
h=figure(2);
plot(a(:,1),1./a(:,2),'ko-','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
axis([0 0.00395,-0.00005 0.001000001])
set(gca,'xtick',0:0.0005:0.0040)
xlabel("\fontsize{27} p");
ylabel("\fontsize{27}\nu");
xlim([0.001 0.0025])
% plot([0.7 0.7],[-1*10^(-4) 1.2*10^(-3)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[-1*10^(-4) 1.2*10^(-3)],'b--','LineWidth',1)
% axis([-0.2 3.2,-1*10^(-4) 1.2*10^(-3)])
