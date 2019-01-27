% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1 - Question 2

function question2()
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

    % Function
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
end
