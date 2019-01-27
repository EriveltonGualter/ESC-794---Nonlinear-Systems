% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1

clear all
close all; clear all

%
question2
question3
question4
question5
%%
% Question 1

% Question 2

A = 1.5;
V1 = 1;
V2 = 1;
L = 1;  % H
C = 1;  % F
R = .1; % Ohm

param = [A V1 V2 L C R];
fun = @(t,x) derStableOscilator(t,x, param);

x0 = [0; 25];
[t, X] =ode45(fun,[0 1000],x0);  

figure('Name','Question 2 - Phase Plane XY','NumberTitle','off');
hold on;
axis equal
plot(X(:,1), X(:,2), 'c'); 
plot(X(1,1), X(1,2), '*r');
plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);

title('Phase Plane XY');
ylabel('x_2'); xlabel('x_1');
legend('','Initial State','Final State');

axis([-30 30 -30 30])

print('question2a','-depsc')

%% Question 3

figure('Name','Question 3 - Phase Portrait','NumberTitle','off');
hold on;
axis equal

xl = 1; yl = 1;
for x = -xl:0.1:xl
    for y = -yl:0.1:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [t, X] =ode45(@derEquation3,[0 20],x0);  
            plot(X(:,1), X(:,2), 'c'); 
            plot(X(1,1), X(1,2), '*r');
            plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);
        end
    end
end

title('Phase Portrait');
ylabel('x_2'); xlabel('x_1');
legend('','Initial State','Final State');
print('question3','-depsc')

J1 = jacob([0 0]);
J2 = jacob([-1 2]);
J3 = jacob([1 2]);

[V1,D1] = eig(J1);
[V2,D2] = eig(J2);
[V3,D3] = eig(J3);

jordan(D1)
jordan(D2)
jordan(D3)

figure; hold on;

line([0 0], [-2.5 2.5], 'color', 'k');
line([-1.5 2.5], [0 0], 'color', 'k');
p1 = plot(real(D1(1,1)), imag(D1(1,1)), '*r', real(D1(2,2)), imag(D1(2,2)), '*r');
p2 = plot(real(D2(1,1)), imag(D2(1,1)), 'sg', real(D2(2,2)), imag(D2(2,2)), 'sg');
p3 = plot(real(D3(1,1)), imag(D3(1,1)), '*b', real(D3(2,2)), imag(D3(2,2)), '*b');

legend([p1(1) p2(1) p3(1)],'eig 1','eig 2','eig 3');
title('Eigenvalues');
ylabel('Imag'); xlabel('Real');
print('question3b','-depsc')

%% Question 4

A = [0 1; -10 -10];
B = [0; 0];
C = [1 0; 0 1];
D = 0;

sys = ss(A,B,C,D);

figure('Name','Question 4 - Phase Portrait','NumberTitle','off');
hold on;
axis equal

xl = 3; yl = 3;
for x = -xl:0.5:xl
    for y = -yl:0.5:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [~,T,X] = initial(sys, x0);
            plot(X(:,1), X(:,2), 'c'); 
            plot(X(1,1), X(1,2), '*r');
            plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);
        end
    end
end

title('Phase Portrait');
ylabel('x_2'); xlabel('x_1');
legend('','Initial State','Final State');
axis(1.5*[-xl xl -yl yl])

print('question4','-depsc')


%% Question 5
figure('Name','Question 5 - Phase Portrait','NumberTitle','off');
hold on;
axis equal

xl = 4; yl = 4;
for x = -xl:1:xl
    for y = -yl:1:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [t, X] =ode45(@derEquation5,[0 20],x0);  
            plot(X(:,1), X(:,2), 'c'); 
            plot(X(1,1), X(1,2), '*r');
            plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);
        end
    end
end

xye = [0 7/3 -7/3; 0 0 0];
for i=1:length(xye)
    [t, X] =ode45(@derEquation5,[0 20],xye(:,i));  
    plot(X(:,1), X(:,2), 'r'); 
end

axis(1.5*[-xl xl -yl yl])

title('Phase Plane XY');
ylabel('x_2'); xlabel('x_1');
legend('','Initial State','Final State');
print('question5','-depsc')


%% Functions
function xDer = derStableOscilator(t, x, param)

    A = param(1);
    V1 = param(2);
    V2 = param(3);
    L = param(4);
    C = param(5);
    R = param(6);
    
    x1 = x(1);
    x2 = x(2);
    
    xd1 = inv(C)*x2;
    xd2 = inv(L)*(V1*atan(A*R*x2/V2)-x1-R*x2) ;
    
    xDer = [xd1; xd2];
end

function xDer = derEquation3 (t, x)
    x1 = x(1);
    x2 = x(2);
    
    xd1 = 2*x1 - x1*x2;
    xd2 = 2*x1*x1 - x2;
    
    xDer = [xd1; xd2];
end

function J = jacob(xe)
    x1e = xe(1); 
    x2e = xe(2); 
    J = [2-x2e, -x1e; 4*x1e, -1];
end

function xDer = derEquation5 (t, x)
    x1 = x(1);
    x2 = x(2);
    
    xd1 = x2;
    xd2 = x1 - 2*atan(x1+x2);
    
    xDer = [xd1; xd2];
end
