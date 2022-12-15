function [dx1,dy1] = rk4(d2x, d2y, dx0, dy0, dt)
    
k1 = d2x(dx0,dy0);
l1 = d2y(dx0,dy0);
k2 = d2x(dx0+.5*dt*k1,dy0+.5*dt*l1);
l2 = d2y(dx0+.5*dt*k1,dy0+.5*dt*l1);
k3 = d2x(dx0+.5*dt*k2,dy0+.5*dt*l2);
l3 = d2y(dx0+.5*dt*k2,dy0+.5*dt*l2);
k4 = d2x(dx0+dt*k3,dy0+dt*l3);
l4 = d2y(dx0+dt*k3,dy0+dt*l3);

dx1 = dx0 + dt*( k1 + 2*k2 + 2*k3 + k4) / 6;
dy1 = dy0 + dt*( l1 + 2*l2 + 2*l3 + l4) / 6;
    
end