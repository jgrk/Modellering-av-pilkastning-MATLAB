clc;clear;close all;

%Uppgift b)

konst.Kx = .001; 
konst.Ky = .01;
konst.h = 1.85;
konst.bulsy = 1.83;
konst.m = 0.026;
konst.g = 9.82;
konst.d = 2.37;
konst.tol = 10^-7;
konst.phi = 2;


V0 = sekmet(@f, 17, 18, konst);

disp("Hastighet som krävs : "+V0)





%Funktioner---------------------------------------------------------

function r = sekmet(f,x1,x0,konst)
t=1;

while abs(t) > konst.tol
    t = f(x1, konst) * ( x1 - x0 ) / ( f(x1, konst) - f(x0, konst) );
    x2 = x1 - t;
    x0 = x1;
    x1 = x2;
    disp([x2 t])

end
r = x2;

end

%-------------------------------------------------------------------

function trff = f(V0, konst)
d2x=@(dx,dy) (- ( konst.Kx / konst.m )* dx* sqrt( dx^2 + dy^2 ) );
d2y=@(dx,dy) (- konst.g-( konst.Ky / konst.m )* dy* sqrt( dx^2 + dy^2 ) );

dt= 10^-6;

clear x y dx dy t

t(1) = 0;
x(1) = 0;
y(1) = konst.h;
dx(1) = V0* cos ( konst.phi* 2* pi / 360 );
dy(1) = V0* sin ( konst.phi* 2* pi / 360 );


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


