% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1

clear all
close all; clear all

% Question 1

% Question 2

% Question 3

A = [0 1; -10 -10];
B = [0; 0];
C = [1 0; 0 1];
D = 0;

sys = ss(A,B,C,D);

figure;
hold on;
axis equal

xl = 1; yl = 1;
for x = -xl:0.1:xl
    for y = -yl:0.1:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [~,T,X] = initial(sys, x0);
            plot(X(:,1), X(:,2), 'b')        
        end
    end
end
axis([-xl xl -yl yl])

%% Question 4

figure;
hold on;
axis equal

xl = 1; yl = 1;
for x = -xl:0.1:xl
    for y = -yl:0.1:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [t, X] =ode45(@derEquation4,[0 20],x0);  
            plot(X(:,1), X(:,2), 'c'); 
            plot(X(1,1), X(1,2), '*r');
            plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);
        end
    end
end

%% Question 5

figure;
hold on;
axis equal

xl = 4; yl = 4;
for x = -xl:1:xl
    for y = -yl:1:yl
        if abs(x) == xl || abs(y) == yl
            x0 = [x; y];
            [t, X] =ode45(@derEquation5,[0 20],x0);  
            plot(X(:,1), X(:,2), 'c'); 
%             plot(X(1,1), X(1,2), '*r');
%             plot(X(end,1), X(end,2), 'ok', 'LineWidth', 2);
        end
    end
end
axis([-xl xl -yl yl])


function xDer = derEquation4 (t, x)
    x1 = x(1);
    x2 = x(2);
    
    xd1 = 2*x1 - x1*x2;
    xd2 = 2*x1*x1 - x2;
    
    xDer = [xd1; xd2];
end

function xDer = derEquation5 (t, x)
    x1 = x(1);
    x2 = x(2);
    
    xd1 = x2;
    xd2 = x1 - 2*atan(x1+x2);
    
    xDer = [xd1; xd2];
end
