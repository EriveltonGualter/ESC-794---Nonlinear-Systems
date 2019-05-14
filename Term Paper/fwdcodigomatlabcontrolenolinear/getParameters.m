function param = getparameters()
    % Parameters

    r = 0.254;  % m
    l = 0.317;  % m
    mw = 29;    % kg
    mb = 310.6; % kg
    Ib = 65;    % Kg.m^2
    Iw = 0.6;   % Kg.m^2
    Dw = 4;     % N.s/m
    Db = 0.1;   % N.s/m
    g = 9.81;   % m/s^2

    param = [r l mw mb Ib Iw Dw Db g];
end