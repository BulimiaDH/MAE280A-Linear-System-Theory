%Problem3
%Assign values to A,B,C,D using values provided by the thesis

A=[-13.692 13.692 128.381;21.023 -21.023 -83.514;1 0 0];
B=[-74.101;113.775;0];
C=[1 0 0;0 1 0];
D=[0;0];
sys=ss(A,B,C,D);
eig(A);

%Using place function to design a pole placement controller for the system
K=place(A,B,[-26.5+i;-26.5-i;26.5])
sys_ctrl=ss(A-B*K,B,C,D);
t=linspace(0,2);
u=zeros(size(t));
lsim(sys,u,t)




