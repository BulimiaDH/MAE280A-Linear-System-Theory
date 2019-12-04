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
I_w=m_w*r^2/2+G_r^2*I_m;

a_val=2*I_w+(m_b+2*m_w)*r^2;
b_val=m_b*r*l;
c_val=I_b+m_b*l^2;
d_val=m_b*g*l;
e_val=2*G_r*s_bar/V_max;
j_val=2*G_r^2*k;

A1_theta=[0 1;subs(subs(theta_A,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(1))]
B1_theta=[0;subs(subs(theta_B,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(1))]
C1_theta=[1 0;0 1]
D1_theta=[0;0];
[num1,den1]=ss2tf(double(A1_theta),double(B1_theta),C1_theta,D1_theta)

A2_theta=[0 1;subs(subs(theta_A,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(2))]
B2_theta=[0;subs(subs(theta_B,[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]),theta,theta_eq(2))]
C2_theta=[1 0;0 1];
D2_theta=[0;0];
[num2,den2]=ss2tf(double(A2_theta),double(B2_theta),C2_theta,D2_theta)

f1=figure('Name','Transfer Function from V to theta_dot about equilibrium 1');
bode(tf(num1(1,3),den1));
f2=figure('Name','Transfer Function from V to theta_dot about equilibrium 2');
bode(tf(num2(1,3),den2));

%Part two
a_val2=2*I_w;
syms phi_d phi_dd;
phi_dd=-j/a*phi_d+e/a*V;
A1_phi=subs(jacobian(phi_dd,phi_d),[a,e,j],[a_val2,e_val,j_val]);
B1_phi=subs(jacobian(phi_dd,V),[a,e,j],[a_val2,e_val,j_val]);
C1_phi=1;
D1_phi=0;
[num3,den3]=ss2tf(double(A1_phi),double(B1_phi),C1_phi,D1_phi)
f3=figure('Name','Transfer Function from V to phi_dot');
bode(tf(num3,den3));



