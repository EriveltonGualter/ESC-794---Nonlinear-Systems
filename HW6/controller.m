function u = controller(X, Kp, Kd)
    ydddot = 0;
    yddot = 0;
    yd = 1;
    
    x1 = X(1);
    x2 = X(2);
    
    e = (X(1) - yd);
    edot = (X(2) - yddot);
    
    v = ydddot - Kd*edot - Kp*e;
    u = v/cos(x2) - x1^4*cos(x2);
end