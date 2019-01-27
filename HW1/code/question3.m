% Erivelton Gualter
% 01/17/2019
%
% Nonlinear System - Homework 1 - Question 3

function question3()

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

    j1 = jordan(J1);
    j2 = jordan(J2);
    j3 = jordan(J3);
    
    D1 = eig(J1);
    D2 = eig(J2);
    D3 = eig(J3);

    figure; hold on;

    line([0 0], [-2.5 2.5], 'color', 'k');
    line([-1.5 2.5], [0 0], 'color', 'k');
    p1 = plot(real(D1(1)), imag(D1(1)), '*r', real(D1(2)), imag(D1(2)), '*r');
    p2 = plot(real(D2(1)), imag(D2(1)), 'sg', real(D2(2)), imag(D2(2)), 'sg');
    p3 = plot(real(D3(1)), imag(D3(1)), '*b', real(D3(2)), imag(D3(2)), '*b');

    legend([p1(1) p2(1) p3(1)],'eig 1','eig 2','eig 3');
    title('Eigenvalues');
    ylabel('Imag'); xlabel('Real');
    print('question3b','-depsc')
    
    figure;
    subplot(121); plotZphase(j1); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    subplot(122); plotZphase(J1); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    
    figure;
    subplot(121); plotZphase(j2); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    subplot(122); plotZphase(J2); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    
    figure;
    subplot(121); plotZphase(j3); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    subplot(122); plotZphase(J3); axis(5*[-1 1 -1 1]); ylabel('z_2'); xlabel('z_1');
    
    
    % Functions
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

    function plotZphase(A)
        sys = ss(A,[0;0],eye(2),0);

        hold on;
        axis equal

        xl = 3; yl = 3;
        for x = -xl:1:xl
            for y = -yl:1:yl
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

    end

end