clc;
clear;

dif = [40; 40];
ms = 200;
xmin = [0; 0];
xmax = [150; 150];
x = linspace(xmin(1), xmax(1), ms);
y = linspace(xmin(2), xmax(2), ms);
[X, Y] = meshgrid(x, y);

file_path = 'pp12_C11_0.030000Xm80.0.txt';
sample=sprintf(file_path);
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
U_max =16.3;%U_max =8.6;
U=U.*(U<U_max)+U_max.*(U>U_max);
pcolor(X, Y, U);
shading interp;
colormap(jet(256));
hold on;
sum(sum(p))
h = colorbar;
h.FontSize = 25;
h.LineWidth = 2;
axis([0 130 ,0 130])

% Calculate the flux line and gradient force line
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
xlabel('\fontsize{25} N1');
ylabel('\fontsize{25} N2');
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])

% Draw normalized flux lines
mg = 1:10:200;
ng = mg;
quiver(X(mg,ng),Y(mg,ng),FFx(mg,ng),FFy(mg,ng),0.45,'color','k','LineWidth',1);
hold on
quiver(X(mg,ng),Y(mg,ng),JJx(mg,ng),JJy(mg,ng),0.5,'color','W','LineWidth',1);

