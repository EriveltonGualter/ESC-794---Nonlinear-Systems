ho = 0.1;   % Initial height, m
hf = 0.6;   % Final height, m

f = @(t,h) tank(t, h, hf, 1);
[ts,hs] = ode45(f,[0 10],ho);
plot(ts, hs);  hold on;

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