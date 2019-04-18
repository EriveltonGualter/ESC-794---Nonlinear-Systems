function dxdt = fcnODE(X,u)
    x1 = X(1);
    x2 = X(2);
    
    xd1 = sin(x2);
    xd2 = x1^4*cos(x2)+u;
    
    dxdt = [xd1; xd2];
end