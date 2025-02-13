clc
clear
a1=0.125:0.0025:0.14;



% a1=0.1:0.1:0.1;
% x0=[64;7;38];  %猜测的初始大概零点值
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-20);%设置求解非线性方程组的精确程度
n_varpara=length(a1);

% 打开一个文件用于保存结果
resultFile = fopen('results.txt', 'w');

for k=1:n_varpara             
    k
    i = 1; 
    x0=[46.284067803441580;41.155887852864076;25.075468667632137]; 
%     x0=xf1;
func=@(x)force_flower(x,a1(k)); 
 xf1=fsolve(func,x0,options)  

 x0=xf1;
  % 将结果保存到文本文件
fprintf(resultFile, ' %.5f  ', a1(k));
fprintf(resultFile, '%.5f ', xf1);
fprintf(resultFile, '\n');

fixfoint(k,:)=xf1;

    % 更新 x0 为 xf1，以便下一次迭代使用
%     x0 = xf1;
end

% 关闭文件
fclose(resultFile);

function f=force_flower(x,a1)
f=zeros(length(x),1);   

N1=x(1);
N2=x(2);
N3=x(3);
d=1;                                              
r1=1/d ;                                           
r2=1/d ;                                           
r3=1/d ;                                           
S1=6  ;                                            
S2=10   ;                                          
S3=14     ;                                        
m1=0.25/d   ;                                      
m2=0.25/d    ;                                     
m3=0.25/d    ;                                     
D=0.25/d  ;                                        
Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1]  ;      % Fig.2A
% Cji=[0.07,0.04,0.04;0.08,0.10,0.08;0.10,0.10,0.14]%Fig.2A
Cji=[0.04,0.07,0.04;0.08,0.08,0.10;0.14,0.10,a1];%Fig.2B
% Cji=[0.04,0.04,0.07;0.10,0.08,0.08;0.10,0.14,0.10]%Fig.2C
%半饱和常数Kij矩阵
% Kji=[1,0.9,0.3;0.3,1,0.9;0.9,0.3,1];  %Fig.3A
% Kji=[1,0.5,0.3;0.3,1,0.5;0.5,0.3,1];  %Fig.3B
% Kji=[1,0.4,0.3;0.3,1,0.4;0.4,0.3,1];  %Fig.3C
%资源含量Cij矩阵
% Cji=[0.04,0.07,0.04;0.08,0.08,0.1;0.14,0.10,0.10];%Fig.3A Fig.3B Fig.3C


R1= S1-Cji(1, 1)*N1-Cji(1, 2)*N2-Cji(1, 3)*N3;    
R2= S2-Cji(2, 1)*N1-Cji(2, 2)*N2-Cji(2, 3)*N3;    
R3= S3-Cji(3, 1)*N1-Cji(3, 2)*N2-Cji(3, 3)*N3;    
P11=r1*R1/(Kji(1,1)+R1) ;                          
P12=r2*R1/(Kji(1,2)+R1) ;                          
P13=r3*R1/(Kji(1,3)+R1)  ;                         
P21=r1*R2/(Kji(2,1)+R2)  ;                         
P22=r2*R2/(Kji(2,2)+R2)   ;                        
P23=r3*R2/(Kji(2,3)+R2)   ;                        
P31=r1*R3/(Kji(3,1)+R3)    ;                       
P32=r2*R3/(Kji(3,2)+R3)    ;                       
P33=r3*R3/(Kji(3,3)+R3)     ;                      
f(1)=N1*(min(min(P11,P21),P31)-m1)  ;             
f(2)=N2*(min(min(P12,P22),P32)-m2)  ;               
f(3)=N3*(min(min(P13,P23),P33)-m3)  ;

end