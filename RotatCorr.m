function newsatpos = RotatCorr(satpos, dTflightSeconds,sat_sys)

We1 = 7.2921151467e-5;      % Earth rotation rate for GPS, QZSS and Galileo  (rad/sec)
We2 = 7.292115e-5;          % Earth rotation rate for Beidou

if sat_sys == "GPS" || sat_sys == "QZSS" || sat_sys == "Galileo"
    we = We1;
elseif sat_sys == "GLONASS" || sat_sys == "BeiDou"
    we = We2;
end


theta = we * dTflightSeconds;

R3 = [ cos(theta)    sin(theta)   0;
      -sin(theta)    cos(theta)   0;
       0                0         1];
   
newsatpos = (R3*satpos(:))';