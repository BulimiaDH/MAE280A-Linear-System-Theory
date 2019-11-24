syms V theta theta_d theta_dd;
theta_dd=1/c*(d*sin(theta)-j*theta_d-e*V);
theta_A=jacobian(theta_dd,[theta,theta_d]);
theta_B=jacobian(theta_dd,V);

theta_eq=[0 pi];
%Substitute parameter values from Appendix A
I_m=3.6e-8;
m_b=263e-3;
m_w=27e-3;
r=34e-3;
l=36e-3;
I_b=4e-4;
g=9.8;
G_r=35.57;
s_bar=0.003;
k=1.25e-6;
V_max=7.4;

a_val=2*I_m+(m_b+2*m_w)*r^2;
b_val=m_b*r*l;
c_val=I_b+m_b*l^2;
d_val=m_b*g*l;
e_val=2*G_r*s_bar/V_max;
j_val=2*G_r^2*k;

A1=[0 1;subs(subs(theta_A,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(1))]
B1=[0;subs(subs(theta_B,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(1))]
C1=[1 0;0 1]
D1=[0;0];
[num1,den1]=ss2tf(double(A1),double(B1),C1,D1)

A2=[0 1;subs(subs(theta_A,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(2))]
B2=[0;subs(subs(theta_B,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(2))]
C2=[1 0;0 1]
D2=[0;0];
[num2,den2]=ss2tf(double(A2),double(B2),C2,D2)

f1=figure('Name','Transfer Function from V to theta_dot about equilibrium 1');
bode(tf(num1(1,3),den1));
f2=figure('Name','Transfer Function from V to theta_dot about equilibrium 2');
bode(tf(num2(1,3),den2));