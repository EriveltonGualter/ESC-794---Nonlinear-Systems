% Erivelton Gualter
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
                plot(X(:,1), X(:,2), 'b'); 
            end
        end
    end
    
    [x1, x2] = meshgrid(-5:0.5:5, -5:.5:5);
    x1dot = A(1,1).*x1 + A(1,2).*x2;
    x2dot = A(2,1).*x1 + A(2,2).*x2;
    quiver(x1,x2,x1dot, x2dot, 'color', 'cyan')

    title('Phase Portrait');
    ylabel('x_2'); xlabel('x_1');
    axis([-xl xl -yl yl])
    box on

    print('question4','-depsc')
end