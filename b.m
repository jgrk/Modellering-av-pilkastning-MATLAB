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

du=@(u) [u(2); 
    (- ( konst.Kx / konst.m )* u(2)* sqrt( u(2)^2 + u(4)^2 ) ); 
    u(4); 
    (- konst.g-( konst.Ky / konst.m )* u(4)* sqrt( u(2)^2 + u(4)^2 ) )];

maxiter = 50;
dt= 0.0063;

clear x y u t

t(1) = 0;
x0 = 0;
y0 = konst.h;
dx0 = konst.V0* cos ( phi* 2* pi / 360 );
dy0 = konst.V0* sin ( phi* 2* pi / 360 );
u(:,1)= [x0; dx0; y0; dy0 ];


while u(1,end) < konst.d
    
    t(end+1) = t(end) + dt;
    k1 = du( u(:,end) );
    k2 = du( u(:,end) + dt*.5*k1 );
    k3 = du( u(:,end) + dt*.5*k2 );
    k4 = du( u(:,end) + dt*k3 );
    u(:,end+1) = u(:,end) + dt*( k1 + 2*k2 + 2*k3 + k4 )/6;


end

%dt2: steglängd som krävs för att x(end) = 2.37

dt2 = ( konst.d - u(1,end-1) ) / u(2,end-1);
t(end) = t(end-1) + dt2;
u(:,end) = u(:,end-1) + du( u(:,end-1) )*dt2;


trff = u(3,end) - konst.bulsy;
   

end


