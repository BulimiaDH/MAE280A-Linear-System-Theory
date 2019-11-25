%Problem3
%Assign values to A,B,C,D using values provided by the thesis

A=[-13.692 13.692 128.381;21.023 -21.023 -83.514;1 0 0];
B=[-74.101;113.775;0];
C=[1 0 0;0 1 0];
D=[0;0];
sys=ss(A,B,C,D);
eig(A);
[sys_num,sys_gen]=ss2tf(A,B,C,D)


%Using place function to design a pole placement controller for the system
K=place(A,B,[-70,-10.8,-200]);
x0=[pi/8,pi/9,-pi/3];
sys_ctrl=ss(A-B*K,B,C,D);
t=linspace(0,10);
u=sin(2*t);
lsim(sys_ctrl,u,t,x0)




