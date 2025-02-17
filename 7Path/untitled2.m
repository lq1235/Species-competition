clc;
clear;
Xm=1.0;
% index1=3;
% index2=4;
dif =[40;40];
ms = 200;
xmin = [0;0];
xmax = [150;150];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

px = load('pp12.txt');
p = reshape(px(:,3),ms,ms);
sum(sum(p))
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;             
PP = eq(Pi,0)+Pi;    
P_eps=min(min(PP));   
P = P_eps*eq(Pi,0)+Pi;
U = -log(P);
U_max =16;
U=U.*(U<U_max)+U_max.*(U>U_max);
surf(X,Y,U)
hold on
shading interp;
colormap([jet(256)]);
% colorbar;
% caxis([-4.5,8]);
view([-44.54,13.855]);
%绘制三维带网格的景观
mesh(X(1:9:200,1:9:200),Y(1:9:200,1:9:200),U(1:9:200,1:9:200)+0.1,'EdgeColor','black','FaceColor','none','LineWidth',2);
%利用函数插值估计二维主路径处U的高度数据值
%从A态到B态
a=load('1-2.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.35,'m','LineWidth',3,'markersize',10)
hold on

% %从B态到A态
a=load('2-1.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.35,'w','LineWidth',3,'markersize',10)
hold on

%从A态到B态
a=load('2-3.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.35,'m','LineWidth',3,'markersize',10)
hold on


% %从B态到A态
a=load('3-2.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.35,'w','LineWidth',3,'markersize',10)
hold on

%从A态到B态
a=load('1-3.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.35,'m','LineWidth',3,'markersize',10)
hold on


% %从B态到A态
a=load('3-1.txt');
x1 = a(:,2);
y1 = a(:,3);
U_path=interp2(X,Y,U,x1,y1);
plot3(x1,y1,U_path+0.4,'w','LineWidth',3,'markersize',10)
hold on

axis([0 115 ,0 115])
% %可逆路径
% a=load('从无病态到疾病态（可逆）.txt');
% x1 = a(:,3);
% y1 = a(:,4);
% U_path=interp2(X,Y,U,x1,y1);
% plot3(x1,y1,U_path+1.2,'w','LineWidth',3,'markersize',10)
% hold on
% %可逆路径
% a=load('从疾病态到无病态（可逆）.txt');
% x1 = a(:,3);
% y1 = a(:,4);
% U_path=interp2(X,Y,U,x1,y1);
% plot3(x1,y1,U_path+1.2,'w','LineWidth',3,'markersize',10)
% hold on
xlabel('\fontsize{25} N_{1}');
ylabel('\fontsize{25} N_{2}');


zlabel('\fontsize{25} U')
% axis([0 0.25 ,0 1])

% set(gca,'xtick',0:1:3)
% set(gca,'ytick',0:1:3)
% set(gca,'ztick',0:5:15)
set(gca,'LineWidth',1,'Fontsize',20)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
view([20.86999999999999,50.09820548974244])
% zlim([-5 15])