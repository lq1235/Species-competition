clc
clear
% seed=1;rng(seed); %使用它就会使得每次运行的随机数一样，不使用的话就不一样，这里不用

dim=3;
D=0.005;

%区域
xmin=[0,0,0];
xmax=[150,150,150];
% randn (1,dim);
    %半饱和常数Kij矩阵
    Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1];  %Fig.2A
    %资源含量Cij矩阵
%     Cji=[0.07,0.04,0.04;0.08,0.10,0.08;0.10,0.10,0.14]%Fig.2A
%     Cji=[0.04,0.04,0.07;0.10,0.08,0.08;0.10,0.14,0.10];%Fig.2C
    Cji=[0.05775,0.07,0.04;0.08,0.08,0.10;0.14,0.10,0.10]%Fig.2B
    %半饱和常Kij矩阵
    % Kji=[1,0.9,0.3;0.3,1,0.9;0.9,0.3,1];  %Fig.3A
    % Kji=[1,0.5,0.3;0.3,1,0.5;0.5,0.3,1];  %Fig.3B
    % Kji=[1,0.4,0.3;0.3,1,0.4;0.4,0.3,1];  %Fig.3C
    %资源含量Cij矩阵
    % Cji=[0.04,0.07,0.04;0.08,0.08,0.10;0.14,0.10,0.1]%Fig.3A Fig.3B Fig.3C
    %参数取值
    d=1;
    ri=1/d;
    m1=0.25/d;m2=0.25/d;m3=0.25/d;
    S1=6;S2=10;S3=14;
    Rjistar=(m1*Kji)./(ri-m1);

        dt=0.1;
        t=0;
        i=1;
        t0=1000;
        steps=t0/dt;
        %初始条件
        N1=1;N2=1;N3=1;
        while t<=t0
            noise = sqrt(2*D)*randn(1,dim);%每次生成的随机数都不一样，所以放在这里
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

            N1= N1+dt*(N1*(min(Pji(:, 1))-m1))+sqrt(dt)*noise(1,1);
            N2= N2+dt*(N2*(min(Pji(:, 2))-m2))+sqrt(dt)*noise(1,2);
            N3= N3+dt*(N3*(min(Pji(:, 3))-m3))+sqrt(dt)*noise(1,3);

            %边界处理
            if N1<xmin(1,1)
                N1=2*xmin(1,1)-N1;
            elseif N1>xmax(1,1)
                N1=2*xmax(1,1)-N1;
            else
                N1=N1;
            end

            if N2<xmin(1,2)
                N2=2*xmin(1,2)-N2;
            elseif N2>xmax(1,2)
                N2=2*xmax(1,2)-N2;
            else
                N2=N2;
            end
            if N3<xmin(1,3)
                N3=2*xmin(1,3)-N3;
            elseif N3>xmax(1,3)
                N3=2*xmax(1,3)-N3;
            else
                N3=N3;
            end

            t=t+dt;
            % x(i,:)=[t,N1,N2,N3,R1,R2,R3];
            x(i,:)=[t,N1,N2,N3];

            i=i+1;
        end
lightBlue = [0.43, 0.33, 0.99];
lightred = [0.98, 0.34, 0.34];

plot(x(:,1),x(:,2), 'Color', lightred,'LineWidth', 1)
hold on
plot(x(:,1),x(:,3), 'Color', lightBlue,'LineWidth', 1)
hold on
plot(x(:,1),x(:,4),'Color', 'green','LineWidth', 1)
hold on

axis([0 300 ,-5 90])
set(gca,'ytick',0:30:90)
set(gca,'xtick',0:50:300)
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'YTickLabelRotation',0);%46是字体的旋转角度
xlabel('Time','Fontsize',25)
ylabel(' Abundance','Fontsize',25)
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])

% legend('N1','N2','N3')
xlabel('Time')
% ylabel('Population abundance','FontSize',15)