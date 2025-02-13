clc
clear
a=xlsread('S3相图.xlsx');
plot(a(:,1),a(:,7)+1.5,'g-o','LineWidth', 1, 'MarkerFaceColor', 'green', 'MarkerSize', 10)
hold on

b=xlsread('S3相图.xlsx');
lightred = [0.98, 0.34, 0.34]; 
lightBlue = [0.43, 0.33, 0.99]; 

plot(b(:,1), b(:,11), 'Color', lightred, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightred, 'MarkerSize', 10);
hold on

c=xlsread('S3相图.xlsx');

plot(c(:,1),c(:,3), 'Color', lightBlue, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightBlue, 'MarkerSize', 10);
hold on

xlabel('\fontsize{27} S_{3}')
ylabel('\fontsize{27} N2')
xlim([0 25])
ylim([-20 150])
set(gca,'xtick',0:5:25)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])