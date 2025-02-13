clc
clear

% Draw the area where C11 is less than or equal to 0.05 and fill it in light green
colo = [0.68, 0.84, 0.95];
rectangle('Position', [-0.032741151324633, -16.439752915470837, 0.05+0.032741151324633, 146.409752915470837], 'FaceColor', colo,'EdgeColor', 'none');
hold on

a=load('data.txt');
lightred = [0.98, 0.34, 0.34];
plot(a(:,1),a(:,2),'Color', lightred, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightred, 'MarkerSize', 10)
hold on

b=load('data1.txt');
plot(b(:,1),b(:,2),'g-o','LineWidth', 1, 'MarkerFaceColor', 'green', 'MarkerSize', 10)
hold on

c=load('data2.txt');
lightBlue = [0.43, 0.33, 0.99];
plot(c(:,1),c(:,2)+2,'Color', lightBlue, 'LineStyle', '-', 'Marker', 'o', 'LineWidth', 1, 'MarkerFaceColor', lightBlue, 'MarkerSize', 10)
hold on

box on;
xlabel('\fontsize{27} C_{11}')
ylabel('\fontsize{27} N1')
% zlabel('N2')
grid on
% set(gca,'xtick',0:0.0005:0.0040)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlim([-0.032741151324633,0.172232848675367])
ylim([-16.439752915470837,129.9702470845292])
legend('S_{1}','S_{3}','S_{2}','FontSize',10)