%10% of a for P
clc 
clear
data = xlsread('Battery Data.xlsx','Battery Outputs','A2:E89660');
time = data(:,1);
SOC = data(:,2)/100;
OCV = data(:,3);
I = data(:,4)/1000;
U = data(:,5)/1000;
T = 1;
n = size(time);
R0(1) = 0.0028;
R1(1) = 0.0035;
C1(1) = 15000;
a1(1) = -(T-2*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
a2(1) = -(R0(1)*T+R1(1)*T+2*R0(1)*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
a3(1) = -(R0(1)*T+R1(1)*T-2*R0(1)*R1(1)*C1(1))/(T+2*R1(1)*C1(1));
lambda = 0.97;
Umod(1) = U(1);
theta = [(1-a1(1))*OCV(1) ; a1(1) ; a2(1) ; a3(1)];
P = abs(diag(0.1*theta)); 
for i = 2:n 
   if I(i)*I(i-1) < 0
    phi = [1 U(i-1) -I(i) I(i-1)];
    Umod(i) = phi*theta;
    e(i) = U(i) - phi*theta;
    K = (P*transpose(phi))/(lambda + phi*P*transpose(phi));
    P = (P-transpose(K)*transpose(phi)*P)/lambda;
    theta = theta + e(i)*K;
    a1(i) = theta(2);
    a2(i) = theta(3);
    a3(i) = theta(4);
    R0(i) = T*(a3(i)-a2(i))/(T+a1(i));
    R1(i) = 2*(a1(i)*(R0(i)^2)-a3(i)*R0(i))/(a3(i)-a2(i)-2*a1(i)*R0(i));
    C1(i) = T*(1+a1(i))/(2*(1-a1(i))*R1(i)); 
   else
    phi = [1 U(i-1) I(i) I(i-1)];
    Umod(i) = phi*theta;
    e(i) = U(i) - phi*theta;
    K = (P*transpose(phi))/(lambda + phi*P*transpose(phi));
    P = (P-K*phi*P)/lambda;
    theta = theta + e(i)*K;
    a1(i) = theta(2);
    a2(i) = theta(3);
    a3(i) = theta(4);
    R0(i) = T*(a3(i)-a2(i))/(T+a1(i));
    R1(i) = 2*(a1(i)*(R0(i)^2)-a3(i)*R0(i))/(a3(i)-a2(i)-2*a1(i)*R0(i));
    C1(i) = T*(1+a1(i))/(2*(1-a1(i))*R1(i)); 
   end
end
plot(time,U)
hold on
plot(time,Umod)
