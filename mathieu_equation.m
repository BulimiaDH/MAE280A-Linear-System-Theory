
function dxdt = mathieu_equation(t,x)
for alpha=-2:0.01:7
    for beta=0:0.01:6
        d=[x(2);-(alpha+beta*cos(t)*x(1))];
        dxdt=d;
        x0=[0;1];
        [t0,y0]=ode45(@mathieu_equation, [0 10],x0);
        plot(t0,y0)
    end
end
return
end





    ?