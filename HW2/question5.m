% Erivelton Gualter
%
% 02/04/2019

clear all
close all

fun = @(t,z) (eqDer(t, z));

x0 = 0*pi/2;
y0 = pi/2;

hold on; box on;
for x0 = -5:1:5
    for y0 = -5:1:5
        [t, y] = ode45(fun, [0 10], [x0 y0]);
        plot(y(:,1), y(:,2), 'b')
    end
end
    
[x1, x2] = meshgrid(-5:.5:5, -5:.5:5);
x1dot = (x1-x2).*(x1.*x1+x2.*x2-1);
x2dot = (x1+x2).*(x1.*x1+x2.*x2-1);
quiver(x1,x2,x1dot, x2dot,'color', 'red')

axis equal
yaxis([-5 5]);


function out = eqDer(t, z)
    x1 = z(1);
    x2 = z(2);
    
    out(1,1) = (x1-x2)*(x1*x1+x2*x2-1);
    out(2,1) = (x1+x2)*(x1*x1+x2*x2-1);
end