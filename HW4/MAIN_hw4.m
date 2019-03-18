% Erivelton Gualter

clear all
close all

a = 1;

figure('Name','Van der Pol system  - Phase Portrait','NumberTitle','off');
hold on; box on;

fun = @(t, x) (eqDer(t, x, a));
for x10 = -8:1:8
    for x20 = -8:1:8
        [t, y] = ode45(fun, [0 20], [x10 x20]);
        plot(y(:,1),y(:,2), 'b')
        axis equal
    end
end

[x1, x2] = meshgrid(-10:1:10, -10:1:10);
x1dot = x2;
x2dot = -(a*(x1.*x1-1).*x2+x1);
quiver(x1,x2,x1dot, x2dot, 'color', 'red')

xlabel('x1'); ylabel('x2');
yaxis([-8 8]);

print('system1','-depsc')

function out = eqDer(t, x, a)
    x1 = x(1); 
    x2 = x(2);  
    
    x1d = x2;
    x2d = -(a*(x1*x1-1)*x2+x1);   
    
    out = [x1d; x2d];
end
    
    