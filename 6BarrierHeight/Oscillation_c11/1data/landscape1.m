clc;
clear;
%Xm=1.0;
% index1=3;
% index2=4;
% dif =[0.05;0.05];
ms = 200;
xmin = [55;0];
xmax = [75;13];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

% a1 =0.8:0.1:0.8;
% num = length(a1);
% file_path = 'pp34_Xm4.0.txt';
% for j = 1:num

% h=figure(j);

% sample=sprintf(file_path,a1(j));
sample = 'pp12_C11_0.056750Xm80.0.txt';
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

% U_max =20;%U_max =8.6;
% U=U.*(U<U_max)+U_max.*(U>U_max);
% U(U==U_max)=NaN;
surf(X,Y,U)
% mesh(X,Y,U);
% pcolor(X,Y,U);
shading interp;
colormap([jet(256)]);

xlabel('\fontsize{25} N1');
ylabel('\fontsize{25} N2');

min(min(U))
max(max(U))

colorbar;
% caxis([-2,U_max]);%caxis([-2,8.6]);
% caxis([min(min(U)),U_max]);
view([-295.29,62.0187]);
% al=min(min(U(88:155,1:52)))
% al1=min(min(U(1:52,88:155)))

zlabel('\fontsize{25} U')
axis([55 70 ,0 11])
% zlim([-5 20])
set(gca,'xtick',55:5:75)
set(gca,'ytick',0:5:13)
% % set(gca,'ztick',0:5:15)
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])

% print(h, '-r600', '-dpdf', ['C15_', num2str(j),'.pdf']);
% print(h, '-r600', '-depsc', ['C15_', num2str(j),'.eps']);

% end