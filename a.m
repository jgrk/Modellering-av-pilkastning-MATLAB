clc;clear;close all;

%Uppgift a) första utkast

konst.Kx = .001; 
konst.Ky = .01;
konst.h = 1.85;
konst.bulsy = 1.83;
konst.m = 0.026;
konst.V0 = 13;
konst.g = 9.82;
konst.phi = 5;
konst.d = 2.37;
konst.tol = 10^-5;

d2x=@(dx,dy) (- ( konst.Kx / konst.m )* dx* sqrt( dx^2 + dy^2 ) );
d2y=@(dx,dy) (- konst.g-( konst.Ky / konst.m )* dy* sqrt( dx^2 + dy^2 ) );

maxiter = 50;
dt=0.1;

y_trff(1) = 1;

for iter = 1:maxiter
    
    clear x y dx dy t
    t(1) = 0;
    x(1) = 0;
    y(1) = konst.h;
    dx(1) = konst.V0* cos ( konst.phi* 2* pi / 360 );
    dy(1) = konst.V0* sin ( konst.phi* 2* pi / 360 );


    while x(end) < konst.d
        
        t(end+1) = t(end) + dt;
        x(end+1) = x(end) + dx(end)*dt;
        y(end+1) = y(end) + dy(end)*dt;    
        [dx(end+1),dy(end+1)] = rk4(d2x,d2y,dx(end),dy(end),dt);
    
    end
    
    %dt2: steglängd som krävs för att x(end) = 2.37
    dt2 = ( konst.d - x(end-1) ) / dx(end-1);
    x(end) = x(end-1) + dx(end-1) * dt2;
    t(end) = t(end-1) + dt2;
    y(end) = y(end-1) + dy(end-1) * dt2;
    
    y_trff(2) = y(end);

    if abs( y_trff(2) - y_trff(1) ) < konst.tol
        trff = y_trff(2);
        break
    
    else
        y_trff(1) = y_trff(2);

        dt = dt / 2;

    end

end

disp("Svar: "+trff)
plot(t,x,t,y,t(end),konst.bulsy,"x")
legend({"x(t)","y(t)","Bullseye"},"Location","SouthEast")
xlabel("tid [s]")
ylabel("Sträcka [m]")



%Funktioner-------------------------------------------------

%Eulers metod
%Runge-Kutta 4



%%

