% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1 - Question 4

function question4()
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
end