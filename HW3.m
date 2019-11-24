sympref('FloatingPointOutput',true)
syms a b c d e j
syms x1 x2 x3 x4 x1d x2d x3d x4d
syms y1 y2 u

%Problem1 Linearization
%Part 1 Linearization about equilibirium 1
x1d=-1/(a*c-b^2*cos(x3)^2)*(-a*d*sin(x3)+a*e*u+b*e*cos(x3)*u+a*j*x1+b*j*cos(x3)*x1+b^2*cos(x3)*sin(x3)*(x1)^2-a*j*x2-b*j*cos(x3)*x2);
x1d_A=jacobian(x1d,[x1,x2,x3,x4]);
x1d_B=jacobian(x1d,u);

x2d=1/(a*c-b^2*cos(x3)^2)*((c+b*cos(x3))*j*x1+b*c*sin(x3)*x1^2-(c+b*cos(x3))*j*x2-b*d*cos(x3)*sin(x3)+(c+b*cos(x3))*e*u);
x2d_A=jacobian(x2d,[x1,x2,x3,x4]);
x2d_B=jacobian(x2d,u);
eq=[0 0 0 0];
x_eq1=subs([x1;x2;x3;x4],[x1,x2,x3,u],eq);
lin_A1=subs(x1d_A,[x1,x2,x3,u],eq);
lin_B1=subs(x1d_B,x3,x_eq1(3));

lin_A2=subs(x2d_A,[x1,x2,x3,u],eq);
lin_B2=subs(x2d_B,x3,x_eq1(3));
A=[lin_A1;lin_A2;1 0 0 0;0 1 0 0];
B=[lin_B1;lin_B2;0;0];

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

for m=1:2
    for n=1:4
        A(m,n)=subs(A(m,n),[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]);
        B(n,:)=subs(B(n,:),[a,b,c,d,e,j],[a_val,b_val,c_val,d_val,e_val,j_val]);
    end
end
A,B,C=[1 0 0 0;0 1 0 0], D=[0;0]

Ob=obsv(A,C);
Co=ctrb(A,B);

%It turns out that the system is controllable but not all states are
%observable

