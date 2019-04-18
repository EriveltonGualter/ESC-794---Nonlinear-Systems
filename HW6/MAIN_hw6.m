% Erivelton

clear all
close all
clc

%% Question 2
X0 = [1.2; 0];
Kp = 200;
Kd = 1;

fun = @(t, X) (fcnODE2(t, X, Kp, Kd));
[t,X] = ode45(fun,[0 15], X0);   

subplot(211); hold on; box on;
    plot(t, X(:,1), 'LineWidth', 2); 
    plot(t, ones(size(t)), 'LineWidth', 2);
    legend('Output', 'Reference');
    ylabel('Output');
    title('Question 5');
    
subplot(212); box on;
    plot(t, X(:,2), 'LineWidth', 2);
    xlabel('Time, s');
    
print('fig_question2', '-dpng');

%% Question 5
X0 = [0; 0];

fun = @(t, X) (fcnODE5(t, X));
[t,X] = ode45(fun,[0 5], X0);   
figure;
subplot(211); hold on; box on;
    plot(t, X(:,1), 'LineWidth', 2); 
    plot(t, pi/2*ones(size(t)), 'LineWidth', 2); 
    legend('Output', 'Reference');
    ylabel('Theta');
    title('Question 5');
    
subplot(212); box on;
    plot(t, X(:,2), 'LineWidth', 2); 
    xlabel('Time, s');
    ylabel('Angular Velocity');
  
print('fig_question5', '-dpng');

%% Functions
function dxdt = fcnODE2(t, X, Kp, Kd)
    % Unpack
    x1 = X(1);
    x2 = X(2);
    
    % Desired Trajectory
    ydddot = 0;
    yddot = 0;
    yd = 1;
    
    % Error
    e = (X(1) - yd);
    edot = (X(2) - yddot);
    
    % Controller
    global U
    v = ydddot - Kd*edot - Kp*e;
    u = v/cos(x2) - x1^4*cos(x2);
    U = [U u];
    
    % Plant
    xd1 = sin(x2);
    xd2 = x1^4*cos(x2)+u;
    
    dxdt = [xd1; xd2];
end

function dxdt = fcnODE5(t, X)
    
    % References
    delta = pi/2;
    
    % Parameters
    k = 4;
    a1 = 1;
    m = 0.1;
    l = 1;
    ko = .02;
    go = -9.81;

    % Unpack
    x1 = X(1) - delta;
    x2 = X(2);
    
    % Control
    u = -k*sign(a1*x1+x2);
    
    x1d = x2;
    x2d = -(go/l)*sin(x1+delta) - ko*x2/m + u/(m*l^2);
    
    dxdt = [x1d; x2d];
end