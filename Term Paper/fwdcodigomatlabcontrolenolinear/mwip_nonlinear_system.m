function [x_dot] = mwip_nonlinear_system(t,x,flag,z)

    % Unpack parameters
    param = getParameters();
    
    r = param(1);
    l = param(2);
    mw = param(3);
    mb = param(4);
    Ib = param(5); 
    Iw = param(6);
    Dw = param(7);
    Db = param(8);
    g = param(9);

    global u s tau_disturbance_hat tau_disturbance_real

    % Unpack
    q1 = x(1);
    q2 = x(2);
    q1d = x(3);
    q2d = x(4);

    m11 = (mb + mw)*r^2 + Iw;
    m12 = mb*l*r;
    m22 = mb*l^2 + Ib;
    Gb = mb*g*l;

    M(1,1) = m11;
    M(1,2) = m12*cos(q2);
    M(2,1) = m11 + m12*cos(q2);
    M(2,2) = m22 + m12*cos(q2);

    N(1,1) = -m12*sin(q2)*q2d^2;
    N(2,1) = -Gb*sin(q2) - m12*sin(q2)*q2d^2;

    F(1,1) = (Dw + Db)*q1d - Db*q2d;
    F(2,1) = Dw*q1d;

% modeling error

% m11error = 0.1;
% m12error = 0.2;
% m22error = 0.3;
% Gberror = 0.4;
% 
% m11_hat = m11 - m11*m11error;
% m12_hat = m12 - m12*m12error;
% m22_hat = m22 - m22*m22error;
% Gb_hat = Gb - Gb*Gberror;

Ib = 66; % 66 kg/m2 is the moment of inertia of the body about the ankle joint axis
mb = 76 + 2.01;
mb_hat = 100;
Ib_hat = Ib*mb_hat/mb;
l_hat =0.7;

m11_hat = (mb_hat + mw)*r^2 + Iw;
m12_hat = mb_hat*l_hat*r;
m22_hat = mb_hat*l_hat^2 + Ib_hat;
Gb_hat = mb_hat*g*l_hat;
 

% controller design
gamma = 3;
%gamma = 4;

d1 =1.5;
d2 = 1.5;
c_smc = 3;
c1 = 2000;
c2 = 0;
c3 = 0;
c4 = 2000;
phi_smc = 0.01;
Tf =0.5;

%a0 = -0.2; % initial condition of q2
%a1 = 0; % initial condition of q2d

a0 = q2;
a1 = q2d;

a2 = -3*(a0/(Tf^2)) -2*(a1/Tf);
a3 = 2*(a0/(Tf^3)) + (a1/(Tf^2));

K = gamma + m11_hat*d2 + (m11_hat + m12_hat*cos(q2))*d1;

if(t > Tf)
    v = 0;
    v_dot=0;
    v_double_dot=0;
else
    v = a0 + a1*t + a2*(t^2) + a3*(t^3);
    v_dot = a1 + 2*a2*t + 3*a3*(t^2)/2;
    v_double_dot = 2*a2 + 3*a3*t;
end

s = q2d + c_smc*q2 - v_dot - c_smc*v;

if abs(s)>phi_smc
    sats=sign(s);
else
    sats=s;
end

tau_disturbance_hat(1,1) = z(1,1) + c1*q1d + c2*q2d;
tau_disturbance_hat(2,1) = z(2,1) + c3*q1d + c4*q2d;

detM_hat = m11_hat*m22_hat - (m12_hat*cos(q2))^2;


u = (1 / (m11_hat + m12_hat*cos(q2)) ) * ...
    ( m11_hat*Gb_hat*sin(q2) - ...
(m12_hat^2)*(q2d^2)*sin(q2)*cos(q2) + m11_hat*tau_disturbance_hat(2,1) - ...
(m11_hat + m12_hat*cos(q2))*tau_disturbance_hat(1,1) + ...
detM_hat*( c_smc*q2d - v_double_dot - c_smc*v_dot) + K*sats);

%disturbance=2.5*sin(4.7*t);
%disturbance = 0;
disturbance = (mb)*sin(q2)*g*l*sin(2*t + pi/2);
tau = [u;0];
tau_ext = [disturbance;0];

x_dot = [x(3);x(4); M \ ( - N - F  + tau + tau_ext)];

M_hat(1,1) = m11_hat;
M_hat(1,2) = m12_hat*cos(q2);
M_hat(2,1) = m11_hat + m12_hat*cos(q2);
M_hat(2,2) = m22_hat + m12_hat*cos(q2);

N_hat(1,1) = -m12_hat*sin(q2)*q2d^2;
N_hat(2,1) = -Gb_hat*sin(q2) - m12_hat*sin(q2)*q2d^2;

    tau_disturbance_real = tau_ext - (M -M_hat)*[x_dot(3,1);x_dot(4,1)] - (N-N_hat) -F;

end