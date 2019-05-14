function [z_dot] = mwip_nonlinear_observer(t, z, flag, x)

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

    % % source http://jn.physiology.org/content/93/1/189
    % Ib = 66; % 66 kg/m2 is the moment of inertia of the body about the ankle joint axis
    % mb = 76 + 2.01; %mB = 76 kg is body mass excluding the feet, mF = 2.01 kg is the mass of the feet, 
    % l = 0.87; %dB = 0.87 m is the height of the body center of mass (COM) above the ankle joint axis
    % r = 0.203;
    % mw = 7;
    % Iw = 0.1445;
    % Dw = 4;
    % Db = 0.1;
    % g = 9.81;

    m11 = (mb + mw)*r^2 + Iw;
    m12 = mb*l*r;
    m22 = mb*l^2 + Ib;
    Gb = mb*g*l;

    q = [x(1,1); x(2,1)];
    q_dot = [x(3,1); x(4,1)];

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

    mb_hat = 100;
    Ib_hat = Ib*mb_hat/mb;
    l_hat =0.7;

    m11_hat = (mb_hat + mw)*r^2 + Iw;
    m12_hat = mb_hat*l_hat*r;
    m22_hat = mb_hat*l_hat^2 + Ib_hat;
    Gb_hat = mb_hat*g*l_hat;

    detM_hat = m11_hat*m22_hat - (m12_hat*cos(q(2,1)))^2;

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
    %a0 = -0.2; % initial condition of q(2,1)
    %a1 = 0; % initial condition of q_dot(2,1)

    a0 = q(2,1);
    a1 = q_dot(2,1);

    a2 = -3*(a0/(Tf^2)) -2*(a1/Tf);
    a3 = 2*(a0/(Tf^3)) + (a1/(Tf^2));

    K = gamma + m11_hat*d2 + (m11_hat + m12_hat*cos(q(2,1)))*d1;

    if(t > Tf)
        v = 0;
        v_dot=0;
        v_double_dot=0;
    else
        v = a0 + a1*t + a2*(t^2) + a3*(t^3);
        v_dot = a1 + 2*a2*t + 3*a3*(t^2)/2;
        v_double_dot = 2*a2 + 3*a3*t;
    end

    s = q_dot(2,1) + c_smc*q(2,1) - v_dot - c_smc*v;

    if abs(s)>phi_smc
        sats=sign(s);
    else
        sats=s;
    end


    tau_disturbance_hat(1,1) = z(1,1) + c1*q_dot(1,1) + c2*q_dot(2,1);
    tau_disturbance_hat(2,1) = z(2,1) + c3*q_dot(1,1) + c4*q_dot(2,1);

    u = (1 / (m11_hat + m12_hat*cos(q(2,1))) ) * ...
        ( m11_hat*Gb_hat*sin(q(2,1)) - ...
    (m12_hat^2)*(q_dot(2,1)^2)*sin(q(2,1))*cos(q(2,1)) + m11_hat*tau_disturbance_hat(2,1) - ...
    (m11_hat + m12_hat*cos(q(2,1)))*tau_disturbance_hat(1,1) + ...
    detM_hat*( c_smc*q_dot(2,1) - v_double_dot - c_smc*v_dot) + K*sats);


    D1 =m22_hat + m12_hat*cos(q(2,1));
    D2 = m11_hat + m12_hat*cos(q(2,1));


    z_dot(1,1) = (1/detM_hat) * ( ( c2*D2 - c1*D1)* ...
        (m12_hat*(q_dot(2,1)^2)*sin(q(2,1)) + u + c1*q_dot(1,1) + c2*q_dot(2,1) + z(1,1)) - ...
        (c2*m11_hat - c1*m12_hat*cos(q(2,1))) * ...
        (Gb_hat*sin(q(2,1)) + m12_hat*(q_dot(2,1)^2)*sin(q(2,1)) + c3*q_dot(1,1) + c4*q_dot(2,1) +z(2,1)));

    z_dot(2,1) = (1/detM_hat) * ( ( c4*D2 - c3*D1)* ...
        (m12_hat*(q_dot(2,1)^2)*sin(q(2,1)) + u + c1*q_dot(1,1) + c2*q_dot(2,1) + z(1,1)) - ...
        (c4*m11_hat - c3*m12_hat*cos(q(2,1))) * ...
        (Gb_hat*sin(q(2,1)) + m12_hat*(q_dot(2,1)^2)*sin(q(2,1)) + c3*q_dot(1,1) + c4*q_dot(2,1) +z(2,1)));



end