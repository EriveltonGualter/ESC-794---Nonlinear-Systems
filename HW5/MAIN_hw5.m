% Erivelton Gualter dos Santos

clear all
close all
clc

% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = tf('s');
G = 5/(s*(s+1)*(s+2));
nyquist(G)

figure;
r = 0;
fun = @(t,y) odeY(t, y, r);

for yo= -5:1:5
    sim('hw5.slx', 20);
    plot(y, yd); hold on; axis equal
end

figure;
yo = 1;
sim('hw5.slx', 20);
plot(t, y); 
xlabel('Time, s');
ylabel('y');

% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [0 1; -1 0];
B = [0;1];

if (rank(ctrb(A,B)) == length(A))
    disp('Controllable');
    K = place(A,B, [-1 -2]);
else
    disp('Not Controllable');
end

t = 0:0.01:8;   
u = zeros(size(t));   
sys = ss(A,B,eye(2),0);
sys_cl = ss(A-B*K,zeros(2,1),eye(2),0);

figure;
lsim(sys_cl,u,t,[1 1])  

% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [1 1; 1 2];
B = [1; 1];
C = [5 1];
D = 0.1;

s = tf('s');

TF = C*inv(s*eye(2) - A)*B;

% Question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ho = 0.1;   % Initial height, m
hf = 0.6;   % Final height, m

figure;
f = @(t,h) tank(t, h, hf, 1);
[ts,hs] = ode45(f,[0 8],ho);
plot(ts, hs, 'LineWidth',2);  hold on;

f = @(t,h) tank(t, h, hf, 2)
[ts,hs] = ode45(f,[0 8],ho);
plot(ts, hs, 'LineWidth',2);  hold on;

f = @(t,h) tank(t, h, hf, 4)
[ts,hs] = ode45(f,[0 8],ho);
plot(ts, hs, 'LineWidth',2);  hold on;

f = @(t,h) tank(t, h, hf, 10)
[ts,hs] = ode45(f,[0 8],ho);
plot(ts, hs, 'LineWidth',2);  hold on;

legend('1','2','4','10');


function yDer = odeY(t, Y, r)
    y1 = Y(1);
    y2 = Y(2);
    y3 = Y(3);
    
    if (r-y1) > 0 
        u = 1;
    else
        u = -1;
    end
    if y1==0
        u = 0;
    end
    
    yd1 = y2;
    yd2 = y3;
    yd3 = -2*y1-3*y3+5*u;
    
    yDer = [yd1; yd2; yd3];
end

function hDer = tank(t, h, hd, beta)
    r = 0.5;        % Radius of Tank, m
    ro = 0.1;       % Outlet opening radius
    a = pi*ro^2;    % cross-section area of the outlet pipe
    g = 9.81; 
    
    Ah = -2*r*h+h^2;
    deltah = h - hd;
    u = a*sqrt(2*g*h)-Ah*beta*deltah;
    hDer = (u - a*sqrt(2*g*h))/Ah;
end




















