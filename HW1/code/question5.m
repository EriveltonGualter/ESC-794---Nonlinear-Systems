% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1 - Question 5

function question5()

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


    % Function
    function xDer = derEquation5 (t, x)
        x1 = x(1);
        x2 = x(2);

        xd1 = x2;
        xd2 = x1 - 2*atan(x1+x2);

        xDer = [xd1; xd2];
    end
end