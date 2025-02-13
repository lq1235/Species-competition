clc
clear;
%initial condition
N1=1;N2=1;N3=1;
%Half full and constant Kij matrix
Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1]; %Fig.2A
%Resource content Cji matrix
Cji=[0.07,0.04,0.04;0.08,0.10,0.08;0.10,0.10,0.14]%Fig.2A
% Cji=[0.04,0.07,0.04;0.08,0.08,0.10;0.14,0.10,0.1]%Fig.2B
% Cji=[0.04,0.04,0.07;0.10,0.08,0.08;0.10,0.14,0.10]%Fig.2C
%参数取值
d=1;
ri=1/d;
m1=0.25/d;m2=0.25/d;m3=0.25/d;
S1=6;S2=10;S3=14;
%The demand of species i on resource j
Rjistar=(m1*Kji)./(ri-m1)

dt=0.01;
t=0;
i=1;

while t<=100
    R1= S1-Cji(1, 1)*N1-Cji(1, 2)*N2-Cji(1, 3)*N3;
    R2= S2-Cji(2, 1)*N1-Cji(2, 2)*N2-Cji(2, 3)*N3;
    R3= S3-Cji(3, 1)*N1-Cji(3, 2)*N2-Cji(3, 3)*N3;
    if R1<=0
        R1=0;
    end
    if R2<=0
        R2=0;
    end
    if R3<=0
        R3=0;
    end
    Rj=[R1,R1,R1;R2,R2,R2;R3,R3,R3];
    Pji=ri*Rj./(Kji+Rj);
    N1= N1+dt*(N1*(min(Pji(:, 1))-m1));
    N2= N2+dt*(N2*(min(Pji(:, 2))-m2));
    N3= N3+dt*(N3*(min(Pji(:, 3))-m3));
    t=t+dt;
    x(i,:)=[t,N1,N2,N3,R1,R2,R3];
    i=i+1;
end

lightBlue = [0.43, 0.33, 0.99];
lightred = [0.98, 0.34, 0.34];

figure(1)
plot(x(:,1),x(:,2),'Color', lightred,'LineWidth', 1)
hold on
plot(x(:,1),x(:,3),'Color', lightBlue,'LineWidth', 1)
hold on
plot(x(:,1),x(:,4),'Color', 'green','LineWidth', 1)
hold on
axis([0 100 ,-5 90])
set(gca,'ytick',0:15:90)
set(gca,'xtick',0:20:100)
set(gca,'XTickLabelRotation',0);
set(gca,'YTickLabelRotation',0);
xlabel('Time','Fontsize',25)
ylabel(' Abundance','Fontsize',25)
% set(gca,'LineWidth',1.2,'Fontsize',25)
% set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
legend('N_{1}','N_{2}','N_{3}')

figure(2)
plot(x(:,1),x(:,5),'Color', lightred,'LineWidth', 1)
hold on
plot(x(:,1),x(:,6),'Color', lightBlue,'LineWidth', 1)
hold on
plot(x(:,1),x(:,7),'Color', 'green','LineWidth', 1)
hold on
legend('R_{1}','R_{2}','R_{3}')
xlabel('Time','Fontsize',25)
ylabel(' quantity of resource','Fontsize',25)