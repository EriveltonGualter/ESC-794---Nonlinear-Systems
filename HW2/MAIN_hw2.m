% Erivelton Gualter
%
% 02/04/2019

clear all
close all

a = 4;

figure('Name','Question 1 - Phase Portrait','NumberTitle','off');
hold on; box on;

fun = @(t, x) (eqDer(t, x, a));
for x10 = -3:1:3
    for x20 = -3:1:3
        [t, y] = ode45(fun, [0 10], [x10 x20]);
        plot(y(:,1),y(:,2), 'b')
        axis equal
    end
end

[x1, x2] = meshgrid(-3:0.5:3, -3:.5:3);
x1dot = x2;
x2dot = sign(-a*x1-x2);
quiver(x1,x2,x1dot, x2dot, 'color', 'red')

yaxis([-3 3]);

print('question1','-depsc')

function out = eqDer(t, x, a)
    x1 = x(1);
    x2 = x(2);  
    
    x1d = x2;
    x2d = sign(-x2-a*x1);   
    
    out = [x1d; x2d];
end
    
    