clear;close all; clc;

%Uppgift b)

konst.m = 0.026;
konst.Kx=0.001;
konst.Ky=0.01;
konst.g=9.82;
konst.d=2.37;
konst.bulsy=1.83;
konst.h=1.85;
konst.V0=13;
konst.tol=10^-5;

%Main--------------------------------------------------------------------

gissa1 = [4 5]; gissa2 = [80 83];
svar1 = sekmet(@(phi) f(phi,konst), gissa1, konst);
svar2 = sekmet(@(phi) f(phi, konst), gissa2, konst);
disp("Svar 1: "+svar1);
disp("Svar 2: "+svar2)


%Funktioner--------------------------------------------------------------

%Sekantmetod
function r = sekmet(f,guess,konst)

t=1; x1=guess(1); x0=guess(2);
while abs(t) > konst.tol
    
    t= f(x1) * ( x1 - x0 ) / ( f(x1) - f(x0) );
    x2 = x1 - t;
    x0 = x1;
    x1 = x2;

end

r= x2;

end


%Räknar utt träffpunkt som en funktion av kastvinkeln phi, enl. uppg. a)
function trff = f(phi, konst)

dx0=konst.V0*cos(phi*2*pi/360);
dy0=konst.V0*sin(phi*2*pi/360);


y0=[0 dx0 konst.h dy0];

t_span = [0 2];

opts = odeset("RelTol",konst.tol,"AbsTol",konst.tol,"Events",@(t,y) stopfun(t,y,konst));

[t,val] = ode45(@(t,y) odefun(t,y,konst),t_span,y0,opts);

y=val(:,3);
trff=y(end)-konst.bulsy;

end

%Se uppg. a)

function dxdt = odefun(t,y,konst)

dxdt=zeros(4,1);

dxdt(1) = y(2);
dxdt(2) = -(konst.Kx/konst.m)*y(2)*sqrt(y(2)^2+y(4)^2);
dxdt(3) = y(4);
dxdt(4) = -konst.g-(konst.Ky/konst.m)*y(4)*sqrt(y(2)^2+y(4)^2);

end

%Se uppg. a)

function [value, isterminal, direction] = stopfun(t,y,konst)

value = y(1) >= konst.d;
isterminal = 1;
direction = 0;

end