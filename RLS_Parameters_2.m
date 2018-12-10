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
b0(1) = R0(1);
b1(1) = -R0(1) + T/C1(1) + (T*R0(1))/(C1(1)*R1(1));
a1(1) = T/(C1(1)*R1(1))-1;
lambda = 0.97;
Umod(1) = U(1);
theta = [b0(1); b1(1); a1(1); OCV(1)];
P = abs(diag(0.01*theta));
OCV_mod(1) = OCV(1);
eocv(1) = 0;
for i = 2:n 
    phi = [I(i); I(i-1); (OCV(i-1)-U(i-1)); 1];
    Umod(i) = transpose(theta)*phi;
    L = (P*phi)/(lambda + transpose(phi)*P*phi);
    %P = (P-L*transpose(phi)*P)/lambda;
    P = (1/lambda)*(eye(4)-L*transpose(phi))*P;
    theta = theta + L*(U(i) - transpose(phi)*theta);
    b0(i) = theta(1);
    b1(i) = theta(2);
    a1(i) = theta(3);
    OCV_mod(i) = theta(4);
    eocv(i) = OCV(i) - OCV_mod(i);
    R0(i) = b0(i);
    R1(i) = (b1(i)-a1(i)*b0(i))/(1+a1(i));
    C1(i) = T/(b1(i)-a1(i)*b0(i)); 
end
plot(time,U)
hold on
plot(time,Umod)
