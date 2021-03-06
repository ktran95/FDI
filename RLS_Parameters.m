clc 
clear
load('data.mat');
time = data(:,1);
SOC = data(:,2)/100;
OCV = data(:,3);
I = data(:,4)/1000;
U = data(:,5)/1000;
T = 1;
n = size(time);
R0(1) = 0.0025;
R1(1) = 0.0035;
C1(1) = 20000;
a1(1) = -(T-2*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
a2(1) = -(R0(1)*T+R1(1)*T+2*R0(1)*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
a3(1) = -(R0(1)*T+R1(1)*T-2*R0(1)*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
lambda = 0.97;
Umod(1) = U(1);
theta = [(1-a1(1))*OCV(1) ; a1(1) ; a2(1) ; a3(1)];
P = abs(diag(0.01*theta)); 
OCV_mod(1) = OCV(1);
eocv(1) = 0;
for i = 2:n 
    phi = [1 U(i-1) I(i) I(i-1)];
    Umod(i) = phi*theta;
    e(i) = U(i) - phi*theta;
    K = (P*transpose(phi))/(lambda + phi*P*transpose(phi));
    P = (P-K*phi*P)/lambda;
    theta = theta + e(i)*K;
    a1(i) = theta(2);
    a2(i) = theta(3);
    a3(i) = theta(4);
    OCV_mod(i) = theta(1)/(1-a1(i));
    eocv(i) = OCV(i) - OCV_mod(i);
    R0(i) = T*(a3(i)-a2(i))/(T+a1(i));
    R1(i) = 2*(a1(i)*(R0(i)^2)-a3(i)*R0(i))/(a3(i)-a2(i)-2*a1(i)*R0(i));
    C1(i) = T*(1+a1(i))/(2*(1-a1(i))*R1(i));
end
plot(time,U)
hold on
plot(time,Umod)
