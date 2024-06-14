function [satpos,satclock] = eph2satpos(SOW,PRN,eph,sat_sys,ps)
% Calculate GPS, Galileo, QZSS, Beidou
% Inputs: 
%  SOW      = Second of week
%  PRN      = PRN number
%  eph      = Ephemaris data
%  sat_sys  = GNSS system
%  ps       = Pseudorange
% Outputs:
%  satpos   = 3by1 ECEF Satellite position [m]
%  satclock = Satellite clock bias [sec]

Sat = find(eph(:,2) == PRN);      % Read selected PRN ephemeride 

if isempty(Sat)
    satpos = nan(1,3);
    satclock = nan(1,1);
    return
end

if sat_sys == "BeiDou"
    SOW = SOW - 14;
end

%% Read Ephemeride
% Orbit Parameters
a       = eph(Sat,13).^2;      % Semi-major axis                       (m)              
e       = eph(Sat,11);         % Eccentricity
w0      = eph(Sat,20);         % Argument of perigee                   (rad)
W0      = eph(Sat,16);         % Right ascension of ascending node     (rad)
Wdot    = eph(Sat,21);         % Rate of right ascension               (rad/sec)
i0      = eph(Sat,18);         % Inclination                           (rad)
idot    = eph(Sat,22);         % Rate of inclination                   (rad/sec)
M0      = eph(Sat,9);          % Mean anomaly                          (rad)
delta_n = eph(Sat,8);          % Mean motion rate                      (rad/sec)

% Correction coefficients
Cuc     = eph(Sat,10);         % Argument of perigee (cos)             (rad) 
Cus     = eph(Sat,12);         % Argument of perigee (sine)            (rad)
Crc     = eph(Sat,19);         % Orbit radius        (cos)             (m)
Crs     = eph(Sat,7);          % Orbit radius        (sine)            (m)
Cic     = eph(Sat,15);         % Inclination         (cos)             (rad) 
Cis     = eph(Sat,17);         % Inclination         (sine)            (rad)

% Time
Toe     = eph(Sat,14);         % Time of Ephemeris                     (SOW : sec of GPS week)
GPS_week= eph(Sat,24);         % GPS Week
Ttm     = eph(Sat,30);         % Transmission time of message -604800  (SOW : sec of GPS week)

% Clock
T0_bias = eph(Sat,3);          % Clock Bias                            (sec)
T0_drift= eph(Sat,4);          % Clock Drift                           (sec/sec)
T0_drate= eph(Sat,5);          % Clock Drift rate                      (sec/sec^2)
Tgd     = eph(Sat,28);         % Time Group delay                      (sec)

% Status
SV_health   = eph(Sat,27);     % SV Health
SV_accuracy = eph(Sat,26);     % SV Accuracy
L2_P_flag   = eph(Sat,25);     % L2 P data flag
L2_code     = eph(Sat,23);     % Code on L2 channel
IODC        = eph(Sat,29);     % Issue of Data, Clock
IODE        = eph(Sat,6);      % Issue of Data, Ephemeris

% Constant
We1 = 7.2921151467e-5;         % Earth rotation rate for GPS, QZSS and Galileo  (rad/sec)
We2 = 7.292115e-5;             % Earth rotation rate for Beidou
mu1 = 3.9860050e14;            % earth gravitational constant for GPS and QZSS
mu2 = 3.986004418e14;          % earth gravitational constant for Galileo and Beiduo
c = 2.99792458e8;

if sat_sys == "GPS" || sat_sys == "QZSS"
    We = We1;
    mu = mu1;
elseif sat_sys == "Galileo"
    We = We1;
    mu = mu2;
else
    We = We2;
    mu = mu2;
end


%% Calculate satellite position and Clock bias

if ~any(Toe <= SOW)
    satpos = nan(1,3);
    satclock = nan(1,1);
    return;
end

% Find correct ephemerides
[~,col] = (min(abs(SOW-Toe)));           % Use closest Toe

tr = ps/c;
TOS = SOW - tr;    % Transition time correction

Tk  = TOS - Toe(col);   % Time from ephemeris reference epoch

if Tk > 302400
    Tk = Tk - 604800;
elseif Tk < -302400
    Tk = Tk + 604800;
end


Mk = M0(col) + (sqrt(mu/a(col)^3) + delta_n(col)) * Tk;    % Mean anomaly at Tk 

% Iterative solution for E 
E_old = Mk;
dE = 1;
while (dE > 10e-12)
  EA = E_old + (Mk-E_old+e(col)*sin(E_old))/(1-e(col)*cos(E_old));
  dE = abs(EA-E_old);
  E_old = EA;
end

u = atan2(sqrt(1-e(col)^2)*sin(EA),cos(EA)-e(col)) + w0(col);
r = a(col)*(1-e(col)*cos(EA));
i = i0(col) + idot(col)*Tk;

% Correction for orbital perturbations
del_u_k = Cuc(col)*cos(2*(u)) + Cus(col)*sin(2*(u));  % Argument of latitude correction
del_r_k = Crc(col)*cos(2*(u)) + Crs(col)*sin(2*(u));  % Radius correction
del_i_k = Cic(col)*cos(2*(u)) + Cis(col)*sin(2*(u));  % Inclination correction

u_k = u + del_u_k;     % Corrected argument of latitude
r_k = r + del_r_k;     % Corrected radius
i_k = i + del_i_k;     % Corrected inclination
  
satpos_in = [r_k*cos(u_k) r_k*sin(u_k) 0]';  % Positions in orbital

if sat_sys == "BeiDou" && PRN <= 5
    SIN_5 = -0.0871557427476582;
    COS_5 = 0.9961946980917456;
    W = W0(col) + Wdot(col)*Tk - We*Toe(col);

    satpos_in_g = [satpos_in(1)*cos(W)-satpos_in(2)*cos(i_k)*sin(W),...
        satpos_in(1)*sin(W)+satpos_in(2)*cos(i_k)*cos(W), satpos_in(2)*sin(i_k)];            
    satpos(1) = satpos_in_g(1)*cos(W*Tk)+satpos_in_g(2)*sin(W*Tk)*COS_5+satpos_in_g(3)*sin(W*Tk)*SIN_5;
    satpos(2) = -satpos_in_g(1)*sin(W*Tk)+satpos_in_g(2)*cos(W*Tk)*COS_5+satpos_in_g(3)*cos(W*Tk)*SIN_5;
    satpos(3) = -satpos_in_g(2)*SIN_5+satpos_in_g(3)*COS_5;
else
    W = W0(col) + (Wdot(col)-We)*Tk - (We*Toe(col));
    satpos(1) = satpos_in(1)*cos(W) - satpos_in(2)*sin(W)*cos(i_k);
    satpos(2) = satpos_in(1)*sin(W) + satpos_in(2)*cos(W)*cos(i_k);
    satpos(3) = satpos_in(2)*sin(i_k);
end


dts = T0_bias(col) + T0_drift(col)*(Tk) + T0_drate(col)*(Tk^2);
dts = dts - 2*sqrt(mu*a(col))*e(col)*sin(EA)/c^2;

satclock = dts-Tgd(col);


end