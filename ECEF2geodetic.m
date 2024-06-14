
function geodetic_latlonhgt = ECEF2geodetic(ECEF_XYZ) 

a = 6378137;                           % semi-major axis
f = 1/298.257223563;                   % flattening
b = a*(1-f);                           % semi-minor axis
e = sqrt(f*(2-f));                     % eccentricity
p = sqrt(ECEF_XYZ(1)^2+ECEF_XYZ(2)^2);                     
lat = atan((ECEF_XYZ(3)/p)*(1-e^2)^-1);          
es = 10^-6;  
maxiter = 50; 
iter = 0;   

while(1)
    oldlat = lat;
    iter = iter + 1; 
    N0 = (a^2)/sqrt((a^2*(cos(oldlat))^2)+((b^2)*(sin(oldlat))^2)); 
    hgt = p/cos(oldlat) - N0;  
    lat = atan ((ECEF_XYZ(3)/p)*(1-((e^2)*N0)/(N0+hgt))^-1);
    ea = abs((lat-oldlat)/lat) * 100 ;
    
    if ea <= es || iter >= maxiter
        break
    end
end

lat = lat * 180 / pi;             
lon = atan2(ECEF_XYZ(2),ECEF_XYZ(1)) * 180 / pi;     %
hgt = round(hgt);

geodetic_latlonhgt = [lat,lon,hgt];



