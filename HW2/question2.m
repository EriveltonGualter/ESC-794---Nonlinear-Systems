% Erivelton Gualter
%
% 02/04/2019

clear all
close all

fun = @(t,z) (eqDer(t, z));

x0 = 0*pi/2;
y0 = pi/2;

[t, y] = ode45(fun, [0 10], [x0 y0]);

plot(y(:,1), y(:,2))

function out = eqDer(t, z)
    x = z(1);
    y = z(2);
    
    out(1,1) =  y + x*(x*x + y*y -1)*sin(1/(x*x + y*y -1));
    out(2,1) = -x + y*(x*x + y*y -1)*sin(1/(x*x + y*y -1));
end