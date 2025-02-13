clc
clear;
% half full and constant Kij matrix
Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1];  %Fig.2A
% Resource content Cji matrix
Cji=[0.04,0.04,0.07;0.1,0.08,0.08;0.10,0.14,0.1];%Fig.2C
d=1;
ri=1/d;
m1=0.25/d;m2=0.25/d;m3=0.25/d;
S1=6;S2=10;S3=14;
% Generates random initial conditions in the range 0 to MAX
MAX=200;
for k=1:2000
    k
    N1 = 0;
    N2 = 0;
    N3 = 0;
    while N1 == 0 || N2 == 0 || N3 == 0
        N1 = rand * MAX;
        N2 = rand * MAX;
        N3 = rand * MAX;
    end
    while N1 == 0 || N2 == 0 || N3 == 0
        N1 = rand * MAX;
        N2 = rand * MAX;
        N3 = rand * MAX;
    end
    disp([N1, N2, N3]);
    [N1 N2 N3]
    dt=0.01;
    t=0;
    i=1;
    while t<=500
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
        if N1<=0
            N1=0;
        end
        if N2<=0
            N2=0;
        end
        if N3<=0
            N3=0;
        end
        t=t+dt;
        x(i,:)=[t,N1,N2,N3,R1,R2,R3];
        i=i+1;
    end
    fixpoint(k,:)=[x(end,2),x(end,3),x(end,4)];
end
% Writes fixpoint to the data.txt file
fid = fopen('data.txt', 'w');
if fid == -1
    error('无法创建文件。');
end
for i = 1:size(fixpoint, 1)
    fprintf(fid, '%f %f %f\n', fixpoint(i, :));
end
fclose(fid);
plot(x(:,1),x(:,2),x(:,1),x(:,3),x(:,1),x(:,4))
legend('1','2','3')
xlabel('Time')
ylabel('Population abundance')