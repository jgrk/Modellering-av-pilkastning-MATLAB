clc;clear;close all;

%Uppgift c)

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



function trff = f(V0, konst)

du=@(u) [u(2); 
    (- ( konst.Kx / konst.m )* u(2)* sqrt( u(2)^2 + u(4)^2 ) ); 
    u(4); 
    (- konst.g-( konst.Ky / konst.m )* u(4)* sqrt( u(2)^2 + u(4)^2 ) )];



dt= 10^-5;

clear x y t u

t(1) = 0;
x0 = 0;
y0 = konst.h;
dx0 = V0* cos ( konst.phi* 2* pi / 360 );
dy0 = V0* sin ( konst.phi* 2* pi / 360 );
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


