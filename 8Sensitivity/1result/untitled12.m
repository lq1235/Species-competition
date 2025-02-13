clc
clear
a=xlsread('BH_Change.xlsx');

% 加载数据
%参考点
b0=a(1,6)-a(1,1),c0=a(1,6)-a(1,2);
b=a(2:1:25,6)-a(2:1:25,1);
c=a(2:1:25,6)-a(2:1:25,2);

b_change=b-b0,c_change=c-c0;

% 创建横坐标的类别
X = categorical({'C_{11}','C_{12}','C_{13}','C_{21}','C_{22}','C_{23}','C_{31}','C_{32}','C_{33}','K_{11}','K_{12}','K_{13}','K_{21}','K_{22}','K_{23}','K_{31}','K_{32}','K_{33}','m_{1}','m_{2}','m_{3}','S_{1}','S_{2}','S_{3}'});
X = reordercats(X,{'C_{11}','C_{12}','C_{13}','C_{21}','C_{22}','C_{23}','C_{31}','C_{32}','C_{33}','K_{11}','K_{12}','K_{13}','K_{21}','K_{22}','K_{23}','K_{31}','K_{32}','K_{33}','m_{1}','m_{2}','m_{3}','S_{1}','S_{2}','S_{3}'});

% 计算每组数据的柱状图的宽度
barWidth = 0.35;

% 计算偏移量，使得两组柱状图可以并列显示
offset = barWidth / 2;

% 转换 categorical 类型为 double 数组
X_numeric = double(X);

% 绘制柱状图
bar(X_numeric - offset, b_change, barWidth, 'FaceColor', [0.98, 0.34, 0.34]); % 蓝色柱状图
hold on
bar(X_numeric + offset, c_change, barWidth, 'FaceColor', [0.43, 0.33, 0.99]); % 红色柱状图

hold off

% 设置横坐标和纵坐标标签
set(gca,'ytick',-5:0.5:5.5)
xticks(X_numeric); % 设置横坐标位置
xticklabels(X);    % 设置横坐标标签
    


ylabel("\fontsize{25} \DeltaBarrier")
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度

% 设置图形属性
set(gca,'LineWidth',1.2,'Fontsize',27.4)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlim([-1,26.150000000000002])
ylim([-1.5 1.5])
% 设置横坐标标签的字体大小
ax = gca;
ax.XAxis.FontSize = 10; % 例如将横坐标字体大小设置为20
% 添加图例，并设置图例大小
legend('S_{1}','S_{2}',  'FontSize', 17) % 添加图例，设置字体大小为20