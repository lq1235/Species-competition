clc
clear
% seed=1;rng(seed); %使用它就会使得每次运行的随机数一样，不使用的话就不一样，这里不用

dim=3;
% dt=0.001;
D=60;%D=20;%小噪声下计算效果也不错

%区域
xmin=[0,0,0];
xmax=[150,150,150];
% randn (1,dim);
a1=9.1:0.4:13.9;
b=length(a1)
for g=1:b
    g

    %半饱和常数Kij矩阵
    Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1];  %Fig.2A
    %资源含量Cij矩阵
    % Cji=[0.07,0.04,0.04;0.08,0.10,0.08;0.10,0.10,0.14]%Fig.2A
    Cji=[0.04,0.04,0.07;0.10,0.08,0.08;0.10,0.14,0.10];%Fig.2C
    % Cji=[0.01,0.07,0.04;0.08,0.08,0.10;0.14,0.10,0.10]%Fig.2B
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
    S1=6;S2=a1(1,g),;S3=14;
    Rjistar=(m1*Kji)./(ri-m1);

    %总共模拟的时间序列计算的次数
    for k=1:1000
        k

        dt=0.1;
        t=0;
        i=1;
%         t0=8000;
        t0=800;%之前用的这个计算的比较好，总共计算了1000次，计算100次的效果也是非常不错，计算50次的效果也不错
        steps=t0/dt;
        %初始条件
        N1=0;N2=99.2857;N3=0;
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

           threshold = 50;
%             限制不能跳到另外一个态
           if any(N1>threshold)
               continue;
           end

                        threshold1 = 40;
                        %         限制不能跳到另外一个态
                        if any(N2<threshold1)
                            continue;
                        end
            
                       threshold3 = 50;
%                                 限制不能跳到另外一个态
                       if any(N3>threshold3)
                           continue;
                       end

            t=t+dt;
            % x(i,:)=[t,N1,N2,N3,R1,R2,R3];
            x(i,:)=[t,N1,N2,N3];

            x0(i,:)=[N1,N2,N3];
            x1=x0';

            i=i+1;

        end



      series1=1; series2=2;series3=3;

        %计算自相关
        traj_1=x1(series1,:);%选取一个变量
        [acf, lags]=autocorr(traj_1,'NumLags',steps-1);
        t1=lags(1:1:end)*dt;
        index(g,k)=min(find(acf<0));
        tau_ik(g,k)=trapz(t1(1:index(g,k)-1),acf(1:index(g,k)-1));
        tau=mean(tau_ik,2);
        %计算自相关
        traj_2=x1(series2,:);%选取一个变量
        [acf1, lags]=autocorr(traj_2,'NumLags',steps-1);
        t1=lags(1:1:end)*dt;
        index1(g,k)=min(find(acf1<0));
        tau_ik1(g,k)=trapz(t1(1:index1(g,k)-1),acf1(1:index1(g,k)-1));
        tau1=mean(tau_ik1,2);
        %计算自相关
        traj_3=x1(series3,:);%选取一个变量
        [acf2, lags]=autocorr(traj_3,'NumLags',steps-1);
        t1=lags(1:1:end)*dt;
        index2(g,k)=min(find(acf1<0));
        tau_ik2(g,k)=trapz(t1(1:index2(g,k)-1),acf2(1:index2(g,k)-1));
        tau2=mean(tau_ik2,2);


        %计算互相关
        traj1=x1([series1,series2],1:end);  %选取两个变量
        [xcf,lags1]=crosscorr(traj1(1,:),traj1(2,:),'NumLags',steps-1);
        CXY=xcf(steps:1:2*steps-1);
        CYX=xcf(steps:-1:1);
        deltaC=CXY-CYX;
        t2=lags1(steps:2*steps-1)*dt;
        n=1000;%对于t趋近于0取平均
        DeltaC(g,k)=sqrt(trapz( t2(1:n) ,deltaC(1:n).^2 )/(dt*n));
        DeltaC1=mean(DeltaC,2);

        %计算互相关
        traj2=x1([series1,series3],1:end);  %选取两个变量
        [xcf1,lags1]=crosscorr(traj2(1,:),traj2(2,:),'NumLags',steps-1);
        CXY1=xcf1(steps:1:2*steps-1);
        CYX1=xcf1(steps:-1:1);
        deltaC1=CXY1-CYX1;
        t2=lags1(steps:2*steps-1)*dt;
        n=1000;%对于t趋近于0取平均
        DeltaC11(g,k)=sqrt(trapz( t2(1:n) ,deltaC1(1:n).^2 )/(dt*n));
        DeltaC2=mean(DeltaC11,2);

        %计算互相关
        traj3=x1([series2,series3],1:end);  %选取两个变量
        [xcf2,lags1]=crosscorr(traj3(1,:),traj3(2,:),'NumLags',steps-1);
        CXY2=xcf2(steps:1:2*steps-1);
        CYX2=xcf2(steps:-1:1);
        deltaC2=CXY2-CYX2;
        t2=lags1(steps:2*steps-1)*dt;
        n=1000;%对于t趋近于0取平均
        DeltaC22(g,k)=sqrt(trapz( t2(1:n) ,deltaC2(1:n).^2 )/(dt*n));
        DeltaC3=mean(DeltaC22,2);

    end
        end
    
%保存计算结果
dlmwrite('tau-tau1-tau2-DeltaC1-DeltaC2-DeltaC3.txt', [a1', tau, tau1, tau2, DeltaC1,DeltaC2,DeltaC3], 'delimiter', '\t');


h=figure(1)

plot(a1,tau,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([9.1 9.1],[2.786310365333501,4.026403065333502],'r--','LineWidth',1)
% axis([6.823471828224232,15.681276828224233,2.786310365333501,4.026403065333502])
title('N1')
xlabel('S_{2}','FontSize',27)
ylabel('\tau','FontSize',27);
print(h, '-r600', '-dpdf', ['CSD_N1','.pdf']);

h=figure(2)
plot(a1,tau1,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([9.1 9.1],[4.127885678698867,6.789885678698869],'r--','LineWidth',1)
% axis([8.006187922387722,14.661187922387729,4.127885678698867,6.789885678698869])
title('N2')
xlabel('S_{2}','FontSize',27)
ylabel('\tau','FontSize',27);
print(h, '-r600', '-dpdf', ['CSD_N2','.pdf']);

h=figure(3)
plot(a1,tau2,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([9.1 9.1],[2.028747832329601,6.421047832329603],'r--','LineWidth',1)
% axis([7.66203455072965,14.982534550729648,2.028747832329601,6.421047832329603])
title('N3')
xlabel('S_{2}','FontSize',27)
ylabel('\tau','FontSize',27);
print(h, '-r600', '-dpdf', ['CSD_N3', '.pdf']);

%时间反演对称性破缺
h=figure(4)
plot(a1,DeltaC1,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([9.1 9.1],[67.19262079195958,92.96078079195959],'r--','LineWidth',1)
% axis([7.412132826217083,15.46468282621708,67.19262079195958,92.96078079195959])
title('N1-N2')
xlabel('S_{2}','FontSize',27)
ylabel('\DeltaCC','FontSize',27);
print(h, '-r600', '-dpdf', ['DeltaCC_N1-N2','.pdf']);

%时间反演对称性破缺
h=figure(5)
plot(a1,DeltaC2,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([9.1 9.1],[58.56481817006458,109.8083181700646],'r--','LineWidth',1)
% axis([7.748975920410855,15.069475920410857,58.56481817006458,109.8083181700646])
title('N1-N3')
xlabel('S_{2}','FontSize',27)
ylabel('\DeltaCC','FontSize',27);
print(h, '-r600', '-dpdf', ['DeltaCC_N1-N3','.pdf']);

%时间反演对称性破缺
h=figure(6)
plot(a1,DeltaC3,'k-o','LineWidth',1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
% plot([9.1 9.1],[49.278446156467496,85.8809461564675],'r--','LineWidth',1)
% axis([7.9,15.354225830430588,49.278446156467496,85.8809461564675])
title('N2-N3')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlabel('S_{2}','FontSize',27)
ylabel('\DeltaCC','FontSize',27);
print(h, '-r600', '-dpdf', ['DeltaCC_N2-N3','.pdf']);


figure(7)
plot(x(:,1),x(:,2),x(:,1),x(:,3),x(:,1),x(:,4))
legend('1','2','3')

figure(8)
plot(x(:,2),x(:,3))
xlabel('N_{1}')
ylabel('N_{2}')
xlim([0 150])
ylim([0 150])

%保存数据
filename1 = ['data_S_{2}=', num2str(9.1, '%.3f'), '--S_{2}=', num2str(13.9, '%.3f'), '.mat'];
save(filename1);
