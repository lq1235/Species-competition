clc
clear

a=xlsread('data.xlsx');
plot(a(:,1),a(:,3),'g-o','LineWidth', 1, 'MarkerFaceColor', 'green', 'MarkerSize', 10)
hold on

b=xlsread('data1.xlsx');
lightBlue = [0.43, 0.33, 0.99]; 
plot(b(:,1), b(:,3), 'Color', lightBlue, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightBlue, 'MarkerSize', 10);
hold on

c=xlsread('data2.xlsx');
lightred = [0.98, 0.34, 0.34]; 
plot(c(:,1),c(:,3)+1.5, 'Color', lightred, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightred, 'MarkerSize', 10);
hold on
plot([4.17 4.17],[-100 200],'r--','LineWidth',1)
hold on
plot([4.56 4.56],[-100 200],'r--','LineWidth',1)
hold on

xlabel('\fontsize{27} S_{1}')
ylabel('\fontsize{27} N2')
xlim([0 8])
ylim([-20 120])
% zlabel('N2')
% grid on
% set(gca,'xtick',3.80:0.4:5.01)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% legend('S_{3}','S_{2}','S_{1}','FontSize',10)
