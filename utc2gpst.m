function [GPStow, GPSweek] = utc2gpst(utc_time)

% GPStow = GPS time of week
% GPSweek = GPS week

% GPS epoch start date
gps_start_date = datetime(1980, 1, 6, 0, 0, 0);

% List of all leap second dates
leap_second_dates = [
    datetime(1981, 7, 1); datetime(1982, 7, 1); datetime(1983, 7, 1);
    datetime(1985, 7, 1); datetime(1987, 1, 1); datetime(1989, 1, 1);
    datetime(1990, 1, 1); datetime(1992, 7, 1); datetime(1993, 7, 1);
    datetime(1994, 7, 1); datetime(1995, 1, 1); datetime(1997, 7, 1);
    datetime(1998, 1, 1); datetime(2005, 1, 1); datetime(2008, 1, 1);
    datetime(2012, 7, 1); datetime(2015, 7, 1); datetime(2016, 1, 1)];

% Number of leap seconds up to the given UTC time
GPStow = zeros(size(utc_time));
GPSweek = zeros(size(utc_time));

for i = 1:numel(utc_time)
    % Number of leap seconds up to the given UTC time
    leap_seconds = sum(leap_second_dates <= utc_time(i));

    % Adjust UTC time by adding leap seconds
    gps_time_datetime = utc_time(i) + seconds(leap_seconds);
    % gps_time_datetime = utc_time(i);

    % Calculate GPS week number and seconds of the week
    gps_week = floor(days(gps_time_datetime - gps_start_date) / 7);
    gps_seconds = round(seconds(gps_time_datetime - (gps_start_date + days(gps_week * 7))));

    % Correct if GPS seconds exceed the weekly limit (604800 seconds in a week)
    if gps_seconds < 0
        gps_week = gps_week - 1;
        gps_seconds = gps_seconds + 604800;
    end

    % Store GPS week and GPS seconds
    GPSweek(i) = gps_week;
    GPStow(i) = gps_seconds;
end
end