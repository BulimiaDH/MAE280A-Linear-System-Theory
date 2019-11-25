%Problem 4 Building an observer/state estimator
A=[-13.692 13.692 128.381;21.023 -21.023 -83.514;1 0 0];
B=[-74.101;113.775;0];
C=[1 0 0;0 1 0];
D=[0;0];

L=place(A',C',[-100 -101 -102])';
K=place(A,B,[-70,-10.8,-200]);
At=[A-B*K B*K;zeros(size(A)) A-L*C];
Bt=[B;zeros(size(B))];
Ct=[C zeros(size(C))];
sys_ob=ss(At,Bt,Ct,D);
x0=[pi/8,pi/9,-pi/3];
t=linspace(0,10);
u=sin(2*t);
lsim(sys_ob,u,t,[x0 x0])