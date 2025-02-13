clc;
clear;
%Xm=1.0;
% index1=3;
% index2=4;
% dif =[0.05;0.05];
ms = 200;
xmin = [0;0];
xmax = [150;150];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

% a1 =1:1:57;
% num = length(a1);
% file_path = 'pp12_S3_%0.6f.txt';
% for j = 1:num

% h=figure(j);

% sample=sprintf(file_path,a1(j));
sample = 'pp12.txt';
px = load(sample);
p = reshape(px(:,3),ms,ms);
sum(sum(p))

FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));

Pi = p/z;

%第一种处理方法
PP = eq(Pi,0)+Pi;
P_eps=min(min(PP));
P = P_eps*eq(Pi,0)+Pi;
% %第二种方法
% eps=1.1e-0;
% P=Pi+eps;

U = -log(P);

% U_max =15.5;%U_max =8.6;
% U=U.*(U<U_max)+U_max.*(U>U_max);
% U(U==U_max)=NaN;
surf(X,Y,U)
% mesh(X,Y,U);
% pcolor(X,Y,U);
shading interp;
colormap([jet(256)]);

xlabel('\fontsize{25} N1');
ylabel('\fontsize{25} N2');


% colorbar;
% caxis([-2,U_max]);%caxis([-2,8.6]);
% caxis([min(min(U)),U_max]);
view([-15.562257647565701,49.405782737936313]);
% al=min(min(U(88:155,1:52)))
% al1=min(min(U(1:52,88:155)))

% zlabel('\fontsize{25} U')
% axis([55 73 ,0 13])
% zlim([6 17])
% set(gca,'xtick',0:30:120)
% set(gca,'ytick',0:30:120)
% set(gca,'ztick',5:10:25)
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'YTickLabelRotation',0);%46是字体的旋转角度
% print(h, '-r600', '-dpdf', ['S3_', num2str(j),'.pdf']);
% print(h, '-r600', '-depsc', ['C15_', num2str(j),'.eps']);

% end