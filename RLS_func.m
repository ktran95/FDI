function [R0,R1,C1,K,P,Umod] = RLS_func(R0,R1,C1,Uk1,Ik1,Uk0,Ik0,T,OCV,P)
a1 = -(T-2*R1*C1)/(T+2*R1*C1);
a2 = -(R0*T+R1*T+2*R0*R1*C1)/(T+2*R1*C1);
a3 = -(R0*T+R1*T-2*R0*R1*C1)/(T+2*R1*C1);
lambda = 0.9999;
theta = [(1-a1)*OCV ; a1 ; a2 ; a3];
phi = [1 Uk0 Ik1 Ik0];
Umod = phi*theta;
e = Uk1 - phi*theta;
K = (P*transpose(phi))/(lambda + phi*P*transpose(phi));
P = (P-K*phi*P)/lambda;
theta = theta + e*K;
a1 = theta(2);
a2 = theta(3);
a3 = theta(4);
R0 = T*(a3-a2)/(T+a1);
R1 = 2*(a1*(R0^2)-a3*R0)/(a3-a2-2*a1*R0);
C1 = T*(1+a1)/(2*(1-a1)*R1);
end