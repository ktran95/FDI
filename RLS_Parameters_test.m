clc 
clear
load('data_new.mat');
load('ocvCurve_allChemistry.mat');
time = data(:,1);
I = data(:,4)/1000;
U = data(:,5)/1000;
OCV = data(:,3);
T = 1;
n = size(time);
R0 = 0.002;
R1 = 0.002;
C1 = 20000;
P = diag([0.1038,0.0078,0.000008,0.000008]); 
T = 1;
U_mod = zeros(n);
R0_curve = zeros(n);
R1_curve = zeros(n);
C1_curve = zeros(n);
SOC = zeros(n);
SOC_round = zeros(n);
%OCV = zeros(n);
%OCV(1) = ocv_curve_nmc_new(101,2);
U_mod(1) = U(1);
R0_curve(1) = R0;
R1_curve(1) = R1;
C1_curve(1) = C1;
SOC(1) = 1;
SOC_round(1) = 1;
Q = 36954*3.6;
for i = 2:n 
    SOC(i) = SOC(i-1) - I(i)/Q;
    SOC_round(i) = round(SOC(i),2);
    A = find(ocv_curve_nmc_new(:,1) == SOC_round(i));
    %OCV(i) = ocv_curve_nmc_new(A,2);
    [R0,R1,C1,K,P,Umod] = RLS_func(R0,R1,C1,U(i),I(i),U(i-1),I(i-1),T,OCV(i),P);
    U_mod(i) = Umod;
    R0_curve(i) = R0;
    R1_curve(i) = R1;
    C1_curve(i) = C1;    
end
plot(time,U)
hold on
plot(time,U_mod)
figure
plot(time,R0_curve)
figure
plot(time,R1_curve)
figure
plot(time,C1_curve)
figure
plot(SOC,OCV)