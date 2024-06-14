function [satpos,satclock] = eph2satpos_GLONASS(SOW,PRN,eph,ps,date_nav)

Sat = find(eph(:,2) == PRN);

pos_x = eph(Sat,6)*1000;      % ECEF x position at t(e) [m]
pos_y = eph(Sat,10)*1000;     % ECEF y position at t(e) [m]
pos_z = eph(Sat,14)*1000;     % ECEF z position at t(e) [m]

vel_x = eph(Sat,7)*1000;   % ECEF x velocity at t(e) [m/s]
vel_y = eph(Sat,11)*1000;  % ECEF y velocity at t(e) [m/s]
vel_z = eph(Sat,15)*1000;  % ECEF z velocity at t(e) [m/s]

acc_x = eph(Sat,8)*1000;   % ECEF x accelation at t(e) [m/s^2]
acc_y = eph(Sat,12)*1000;  % ECEF x accelation at t(e) [m/s^2]
acc_z = eph(Sat,16)*1000;  % ECEF x accelation at t(e) [m/s^2]

% Clock
T0_bias = eph(Sat,3);          % Clock Bias [sec]
T0_Fre_bias = eph(Sat,4);      % Frequency Bias

Date = date_nav(Sat);

We = 7.292115e-5;    % Earth angular velocty [rad/s]
mu = 3.9860044e14;   % gravitational constant
J2 = 1.0826257e-3;   % 2nd zonal harmonic of geopot
Re = 6378136;        % Radius of earth [m]
speed_of_light = 2.99792458e8;


if isempty(Sat)
    satpos = nan(1,3);
    satclock = nan(1,1);
    return
end


%% Satellite Position Computation %%

% Find correct ephemerides

Toe = utc2gpst(Date);
SOW = SOW - 18;

[~,col] = min(abs(SOW-Toe));           % Use closest Toe

pos_te = [pos_x(col) pos_y(col) pos_z(col)];
vel_te = [vel_x(col) vel_y(col) vel_z(col)];
acc_te = [acc_x(col) acc_y(col) acc_z(col)];

tr = ps/speed_of_light;
TOS = SOW - tr;
t = TOS - Toe(col);


% t = SOW - Toe(col);

satclock = -T0_bias(col) + T0_Fre_bias(col)*t;

if t < 0
    TT = -60;
else
    TT = 60;
end

% Runge-Kutta
x_init = [pos_te vel_te];
xdot = zeros(1,6);
K = zeros(4,6);      % Runge-Kutta Constant
x = x_init;
sat_pos = x;

while abs(t) > 1e-9
    if abs(t) < 60
        TT = t;
    end

    sat_pos_old = sat_pos;
    x = sat_pos;

    for i = 1:4
        r2 = dot(x(1:3),x(1:3));  % r^2
        r3 = r2*sqrt(r2);         % r^3
        omeg2 = We^2;

        a = 1.5*J2*mu*Re^2/r2/r3;
        b = 5*x(3)^2/r2;
        c = mu/r3-a*(1-b);

        xdot(1) = x(4);
        xdot(2) = x(5);
        xdot(3) = x(6);
        xdot(4) = (c+omeg2)*x(1)+2*We*x(5)+acc_te(1);
        xdot(5) = (c+omeg2)*x(2)-2*We*x(4)+acc_te(2);
        xdot(6) = (c-2*a)*x(3)+acc_te(3);

        K(i,:) = xdot;

        x = sat_pos + K(i,:).*(TT/2);

        if i == 3
            x = sat_pos + K(i,:).*TT;
        end
    end

    for i = 1:6
        sat_pos(i) = sat_pos_old(i) + (K(1,i)+2*K(2,i)+2*K(3,i)+K(4,i))*TT/6;
    end

    t = t - TT;

end

satpos = sat_pos(1:3);

