clc
clear

% dir = 'G:\cell_cycle\2021nature--Spatiotemporal dissection of the cell cycle with single-cell proteogenomics\landscape\multi\0\';
% dir = 'D:\Users\USER\Desktop\Hussman模型\8Coherence\2相干性\备份\';

% dor = 'D=0';
% n_tra = '2';
% filename1= [dir, 'path_time', num2str(n_tra), '.txt'];
% path = importdata(filename1);
a=load("0.1375.txt");
 a(:, 2);
 a(:, 3);
hold on
plot( a(:, 2),  a(:, 3))
% title(dor);
% xlabel('UMAP1');
% ylabel('UMAP2');
% xlim([-5, 12]);
% ylim([0, 12]);
% box on
% set(gca, 'FontName', 'Arial')
% set(gca,'FontSize',20, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1)
% set(gca,'TickDir', 'in', 'TickLength', [0.02 0])
% set(gca, 'LineWidth', 1.5, 'Color', [0 0 0])
% set(gca, 'XColor', [0.00 0.00 0.00])
% set(gca, 'YColor', [0.00 0.00 0.00])
% set(gca, 'ZColor', [0.00 0.00 0.00])
% pbaspect([1 1 1])
set(gca, 'color', 'white');
% set(gca, 'LineWidth', 1.5)
% legend(dor, 'FontSize', 20, 'TextColor', [0 0 0], 'EdgeColor', [0 0 0])
% set(gca, 'color', 'white')
% saveas(figure(1), [dir, 'phasetraj2_', num2str(n_tra),'.fig']); 
% print(figure(1), '-r600', '-dpdf', [dir, 'phasetraj2_', num2str(n_tra),'.pdf']);

n=9999;
x_center = mean(a(5000:n+1,2));
y_center = mean(a(5000:n+1,3));
c0 = x_center + y_center*1i; %center
r1 = a(1:n+1,2:3);



rr = r1(:,1)+r1(:,2)*1i-c0;

ang = angle(rr);
phi = diff(ang);
phi = phi - 2*pi*(phi>pi);
phi = phi + 2*pi*(phi<-pi);

B = sum(abs(phi));
phi0 = phi.*(phi>0);
A = sum(phi0);
coherence = 2*A/B - 1