clc;
clear;

dif = [35; 35];
ms = 200;
xmin = [0; 0];
xmax = [120; 120];
x = linspace(xmin(1), xmax(1), ms);
y = linspace(xmin(2), xmax(2), ms);
[X, Y] = meshgrid(x, y);
a1=1.0:0.1:7.9;
num = length(a1);
file_path = 'pp23_S1_%0.6f.txt';
for j = 1:num
h=figure(j);
sample=sprintf(file_path,a1(j));
px = load(sample);
p = reshape(px(:, 3), ms, ms);
FPx = reshape(px(:, 4), ms, ms);
FPy = reshape(px(:, 5), ms, ms);
z = trapz(y, trapz(x, p));
Pi = p / z;
PP = eq(Pi, 0) + Pi;
P_eps = min(min(PP));
P = P_eps * eq(Pi, 0) + Pi;
U = -log(P);
pcolor(X, Y, U);
shading interp;
colormap(jet(256));
hold on;
sum(sum(p))
dx = x(2) - x(1);
dy = y(2) - y(1);
[GPx, GPy] = gradient(P, dx, dy);
Jx = FPx .* P - dif(1) * GPx;
Jy = FPy .* P - dif(2) * GPy;
E = Jy.^2 + Jx.^2;
JJx = Jx ./ (sqrt(E) + eps);
JJy = Jy ./ (sqrt(E) + eps);

Fx= dif(1)*GPx./P;
Fy=dif(2)*GPy./P;
F=Fx.^2+Fy.^2;
FFx=Fx./(sqrt(F)+eps);
FFy=Fy./(sqrt(F)+eps);
current_aspect_ratio = daspect;
daspect(current_aspect_ratio);
xlabel('\fontsize{25} N2');
ylabel('\fontsize{25} N3');
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
mg = 1:10:200;
ng = mg;
% 
% for i = 1:numel(mg)
%     for j = 1:numel(ng)
%         %计算箭头的起点和终点
%         start_point = [X(mg(i), ng(j)), Y(mg(i), ng(j)), 0];
%         end_point = [X(mg(i), ng(j)) + JJx(mg(i), ng(j)), Y(mg(i), ng(j)) + JJy(mg(i), ng(j)), 0];
%         
%         %确保起点和终点不相同
%         if norm(start_point - end_point) > 0
% 
%            % 绘制箭头
%             arrow3(start_point, end_point, 'w-', 0.8, 2, 0.00001);
%             hold on
%             %绘制线段
%             plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'w-',LineWidth=1);
%             hold on;
% 
%         end
%     end
% end

% 绘制归一化力线
% for i = 1:numel(mg)
%     for j = 1:numel(ng)
%         % 计算力的起点和终点
%         start_point = [X(mg(i), ng(j)), Y(mg(i), ng(j)), 0];
%         end_point = [X(mg(i), ng(j)) + FFx(mg(i), ng(j)), Y(mg(i), ng(j)) + FFy(mg(i), ng(j)), 0];
%         
%         % 确保起点和终点不相同
%         if norm(start_point - end_point) > 0
%             % 绘制力的箭头
%             arrow3(start_point, end_point, 'k-', 0.8, 2, 0.00001);
%             hold on
%              % 绘制线段
%             plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'k-',LineWidth=1);
%             hold on;
%         end
%     end
% end

quiver(X(mg,ng),Y(mg,ng),FFx(mg,ng),FFy(mg,ng),0.45,'color','k','LineWidth',1);
hold on
quiver(X(mg,ng),Y(mg,ng),JJx(mg,ng),JJy(mg,ng),0.5,'color','W','LineWidth',1);
name = 'S1_%0.6f.txt';
sample=sprintf(name,a1(j));
print(h, '-r600', '-dpdf', ['S1_', sample,'.pdf']);
% print(h, '-r600', '-depsc', ['C15_', num2str(j),'.eps']);
end
