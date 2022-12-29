clc;clear;close all;

%Uppgift a)

%Konstanter/Parametrar

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

%Reducerat ekvationssystem
du=@(u) [u(2); 
    (- ( konst.Kx / konst.m )* u(2)* sqrt( u(2)^2 + u(4)^2 ) ); 
    u(4); 
    (- konst.g-( konst.Ky / konst.m )* u(4)* sqrt( u(2)^2 + u(4)^2 ) )];

maxiter = 50;
dt=0.1;

y_trff(1) = 1;

for iter = 1:maxiter
    
    clear x y t u
    
    t(1) = 0;
    x0 = 0;
    y0 = konst.h;
    dx0 = konst.V0* cos ( konst.phi* 2* pi / 360 );
    dy0 = konst.V0* sin ( konst.phi* 2* pi / 360 );
    u(:,1)= [x0; dx0; y0; dy0 ];


    while u(1,end) < konst.d
        
        t(end+1) = t(end) + dt;

        %RK4
        k1 = du( u(:,end) );
        k2 = du( u(:,end) + dt*.5*k1 );
        k3 = du( u(:,end) + dt*.5*k2 );
        k4 = du( u(:,end) + dt*k3 );
        u(:,end+1) = u(:,end) + dt*( k1 + 2*k2 + 2*k3 + k4 )/6;

    
    end
    
    %Interpolation, linjär.
    %Sista steglängden görs med dt2 som är s.a. sista x-elementet är exakt 
    %2.37. 
    
    dt2 = 6*( konst.d - u(1,end-1) ) / (k1(1)+2*k2(1)+2*k3(1)+k4(1));
    t(end) = t(end-1) + dt2;
    u(:,end) = u(:,end-1) + ( k1+2*k2+2*k3+k4 )*dt2/6;

    x = u(1,:);
    y = u(3,:);

    
    y_trff(2) = y(end);

    if abs( y_trff(2) - y_trff(1) ) < konst.tol
        trff = y_trff(2);
        break
    
    else
        y_trff(1) = y_trff(2);

        dt = dt / 2;

    end

end

%Presentation
disp("Svar: "+trff)
plot(t,x,t,y,t(end),konst.bulsy,"x")
legend({"x(t)","y(t)","Bullseye"},"Location","SouthEast")
xlabel("tid [s]")
ylabel("Sträcka [m]")

