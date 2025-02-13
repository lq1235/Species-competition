clc
clear;
%initial condition
N1=1;N2=1;N3=1;
%Half full and constant Kij matrix
Kji=[1,0.6,0.3;0.3,1,0.6;0.6,0.3,1];  %Fig.2A
%Resource content Cji matrix
Cji=[0.07,0.04,0.04;0.08,0.10,0.08;0.10,0.10,0.14]%Fig.2A
% Cji=[0.04,0.07,0.04;0.08,0.08,0.10;0.14,0.10,0.1]%Fig.2B
% Cji=[0.04,0.04,0.07;0.10,0.08,0.08;0.10,0.14,0.10]%Fig.2C
%parameter
d=1;
ri=1/d;
m1=0.25/d;m2=0.25/d;m3=0.25/d;
S1=6;S2=10;S3=14;
%The demand of species i on resource j
Rjistar=(m1*Kji)./(ri-m1)
