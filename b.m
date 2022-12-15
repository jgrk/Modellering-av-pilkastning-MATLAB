clc;clear;close all;

%Uppgift b)

konst.Kx = .001; 
konst.Ky = .01;
konst.h = 1.85;
konst.bulsy = 1.83;
konst.m = 0.026;
konst.V0 = 13;
konst.g = 9.82;
konst.d = 2.37;
konst.tol = 10^-4;


phi1 = sekmet(@f, 4, 5, konst);
phi2 = sekmet(@f, 80, 83, konst);

disp("Vinkel 1 : "+phi1)
disp("Vinkel 2 : "+phi2)




%Funktioner

function r = sekmet(f,x1,x0,konst)
t=1;

while abs(t) > konst.tol
    t = f(x1, konst) * ( x1 - x0 ) / ( f(x1, konst) - f(x0, konst) );
    x2 = x1 - t;
    x0 = x1;
    x1 = x2;

end
r = x2;

end



function trff = f(phi, konst)
d2x=@(dx,dy) (- ( konst.Kx / konst.m )* dx* sqrt( dx^2 + dy^2 ) );
d2y=@(dx,dy) (- konst.g-( konst.Ky / konst.m )* dy* sqrt( dx^2 + dy^2 ) );

maxiter = 50;
dt= 9.765625000000000e-05;



clear x y dx dy t

t(1) = 0;
x(1) = 0;
y(1) = konst.h;
dx(1) = konst.V0* cos ( phi* 2* pi / 360 );
dy(1) = konst.V0* sin ( phi* 2* pi / 360 );


while x(end) < konst.d
    
    t(end+1) = t(end) + dt;
    x(end+1) = x(end) + dx(end)*dt;
    y(end+1) = y(end) + dy(end)*dt;    
    [dx(end+1),dy(end+1)] = rk4(d2x,d2y,dx(end),dy(end),dt);

end

%dt2: steglängd som krävs för att x(end) = 2.37

dt2 = ( konst.d - x(end-1) ) / dx(end-1);
y(end) = y(end-1) + dy(end-1) * dt2;

trff = y(end) - konst.bulsy;
   

end


