clear all
clc
close all

%%
obs = rinexread("Data\YONS00KOR_R_20241650000_01H_01S_MO.rnx");
obs_info = rinexinfo("Data\YONS00KOR_R_20241650000_01H_01S_MO.rnx");
nav = rinexread("Data\YONS00KOR_R_20241650000_01H_01S_MN.rnx");


ref_pos = [-3047506.741 4043980.58 3865243.01];
ref_pos_geodetic = ecef2lla(ref_pos);
speed_of_light = 2.99792458e8 ;


%% GNSS Data
Data_OBS.header = obs_info;
Data_OBS.date = unique(obs.GPS.Time);
epoch = length(Data_OBS.date);
[Data_OBS.GPStow, Data_OBS.GPSweek] = utc2gps(Data_OBS.date);

Data_OBS.pseudorange = nan(epoch,129);

for i = 1:epoch

    % GPS
    idx = obs.GPS.Time == Data_OBS.date(i);
    sat = obs.GPS.SatelliteID(idx);
    GPS_pseudorange_temp = obs.GPS.C1C(idx);
    for j = 1:length(sat)
        Data_OBS.pseudorange(i,sat(j)) = GPS_pseudorange_temp(j);
    end

    % GLONASS
    idx = obs.GLONASS.Time == Data_OBS.date(i);
    sat = obs.GLONASS.SatelliteID(idx);
    GLONASS_pseudorange_temp = obs.GLONASS.C1C(idx);
    for j = 1:length(sat)
        Data_OBS.pseudorange(i,32+sat(j)) = GLONASS_pseudorange_temp(j);
    end

    % BeiDou
    idx = obs.BeiDou.Time == Data_OBS.date(i);
    sat = obs.BeiDou.SatelliteID(idx);
    BeiDou_pseudorange_temp = obs.BeiDou.C2I(idx);
    for j = 1:length(sat)
        Data_OBS.pseudorange(i,61+sat(j)) = BeiDou_pseudorange_temp(j);
    end

    % Galileo
    idx = obs.Galileo.Time == Data_OBS.date(i);
    sat = obs.Galileo.SatelliteID(idx);
    Galileo_pseudorange_temp = obs.Galileo.C1X(idx);
    for j = 1:length(sat)
        Data_OBS.pseudorange(i,90+sat(j)) = Galileo_pseudorange_temp(j);
    end

    % QZSS
    idx = obs.QZSS.Time == Data_OBS.date(i);
    sat = obs.QZSS.SatelliteID(idx);
    QZSS_pseudorange_temp = obs.QZSS.C1C(idx);
    for j = 1:length(sat)
        Data_OBS.pseudorange(i,120+sat(j)) = QZSS_pseudorange_temp(j);
    end
end

%% Satellite Position Data
% PRN(GPS: 1~32, GLONASS: 38~61, BeiDou: 62~90, GALILEO: 91~120, QZSS: 121~129)

for j = 1:5
    if j == 1
        data_obs = obs.GPS;
        data_nav = nav.GPS;
        sat_pos.GPS = nan(epoch,32*3);
        sat_clock.GPS = nan(epoch,32);
        sat_sys = "GPS";
        true_sat_pos.GPS = nan(epoch,32*3);
    elseif j == 2
        data_obs = obs.GLONASS;
        data_nav = nav.GLONASS;
        sat_pos.GLONASS = nan(epoch,24*3);
        sat_clock.GLONASS = nan(epoch,24);
        sat_sys = "GLONASS";
    elseif j == 3
        data_obs = obs.BeiDou;
        data_nav = nav.BeiDou;
        sat_pos.BeiDou = nan(epoch,29*3);
        sat_clock.BeiDou = nan(epoch,29);
        sat_sys = "BeiDou";
    elseif j == 4
        data_obs = obs.Galileo;
        data_nav = nav.Galileo;
        sat_pos.Galileo = nan(epoch,30*3);
        sat_clock.Galileo = nan(epoch,30);
        sat_sys = "Galileo";
    elseif j == 5
        data_obs = obs.QZSS;
        data_nav = nav.QZSS;
        sat_pos.QZSS = nan(epoch,9*3);
        sat_clock.QZSS = nan(epoch,9);
        sat_sys = "QZSS";
    end

    date_nav = data_nav.Time;
    data_nav = timetable2table(data_nav);
    data_nav = table2array(data_nav(:,2:end));
    [s1,s2] = size(data_nav);
    data_nav = [zeros(s1,1) data_nav];


    for i = 1:epoch
        eph_t = find(data_obs.Time == (Data_OBS.date(i)));     % Check time
        PRN   = data_obs.SatelliteID(eph_t);             % Check satellite number (PRN)
        date = Data_OBS.date(i);
        if sat_sys == "BeiDou"
            ps = data_obs.C2I(eph_t);
        elseif sat_sys == "Galileo"
            ps = data_obs.C1X(eph_t);
        else
            ps = data_obs.C1C(eph_t);
        end

        satpos = nan(length(PRN),3);
        satclock = nan(length(PRN),1);
        
        if j == 2
            continue
        end

        for ii = 1:length(PRN)
            if sat_sys == "GLONASS"
                [satpos(ii,1:3),satclock(ii)] = eph2satpos_GLONASS(Data_OBS.GPStow(i),PRN(ii),data_nav,ps(ii),date_nav);
            else
                [satpos(ii,1:3),satclock(ii)] = eph2satpos(Data_OBS.GPStow(i),PRN(ii),data_nav,sat_sys,ps(ii));
            end

            if sat_sys == "GPS"
                sat_pos.GPS(i,3*PRN(ii)-2:3*PRN(ii)) = satpos(ii,1:3);
                sat_clock.GPS(i,PRN(ii)) = satclock(ii);
                % use_nav_true = find(nav.GPS.SatelliteID == PRN(ii));
                % if isempty(use_nav_true)
                %     continue;
                % else
                %     [true_sat_pos.GPS(i,3*PRN(ii)-2:3*PRN(ii)),~,~] = gnssconstellation(Data_OBS.date(i)-seconds(18),nav.GPS(use_nav_true(1),:));
                % end

            elseif sat_sys == "GLONASS"
                sat_pos.GLONASS(i,3*PRN(ii)-2:3*PRN(ii)) = satpos(ii,1:3);
                sat_clock.GLONASS(i,PRN(ii)) = satclock(ii);
            elseif sat_sys == "BeiDou"
                sat_pos.BeiDou(i,3*PRN(ii)-2:3*PRN(ii)) = satpos(ii,1:3);
                sat_clock.BeiDou(i,PRN(ii)) = satclock(ii);
            elseif sat_sys == "Galileo"
                sat_pos.Galileo(i,3*PRN(ii)-2:3*PRN(ii)) = satpos(ii,1:3);
                sat_clock.Galileo(i,PRN(ii)) = satclock(ii);
            elseif sat_sys == "QZSS"
                sat_pos.QZSS(i,3*PRN(ii)-2:3*PRN(ii)) = satpos(ii,1:3);
                sat_clock.QZSS(i,PRN(ii)) = satclock(ii);
            end


        end
    end
end
%%
sat_pos.All = nan(epoch,129*3);
sat_clock.All = nan(epoch,129);

sat_pos.All(:,1*3-2:32*3) = sat_pos.GPS;
sat_pos.All(:,38*3-2:61*3) = sat_pos.GLONASS;
sat_pos.All(:,62*3-2:90*3) = sat_pos.BeiDou(:,1:29*3);
sat_pos.All(:,91*3-2:120*3) = sat_pos.Galileo(:,1:30*3);
sat_pos.All(:,121*3-2:129*3) = sat_pos.QZSS;

sat_clock.All(:,1:32) = sat_clock.GPS;
sat_clock.All(:,38:61) = sat_clock.GLONASS;
sat_clock.All(:,62:90) = sat_clock.BeiDou(:,1:29);
sat_clock.All(:,91:120) = sat_clock.Galileo(:,1:30);
sat_clock.All(:,121:129) = sat_clock.QZSS;

%% Positining
wgs84 = wgs84Ellipsoid('meter');

use_GNSS = "GPS";

if use_GNSS == "GPS"
    pseudorange_value = Data_OBS.pseudorange(:,1:32);
    sat_pos_value = sat_pos.GPS;
    % sat_pos_value = true_sat_pos.GPS;
    % sat_pos_value = Satpos;
    sat_clock_value = sat_clock.GPS;
elseif use_GNSS == "GLONASS"
    pseudorange_value = Data_OBS.pseudorange(:,38:61);
    sat_pos_value = sat_pos.GLONASS;
%     sat_pos_value = Satpos(:,38*3-2:61*3);
    sat_clock_value = sat_clock.GLONASS;
elseif use_GNSS == "BeiDou"
    pseudorange_value = Data_OBS.pseudorange(:,62:90);
    sat_pos_value = sat_pos.BeiDou;
    sat_clock_value = sat_clock.BeiDou;
elseif use_GNSS == "Galileo"
    pseudorange_value = Data_OBS.pseudorange(:,91:120);
    sat_pos_value = sat_pos.Galileo;
    sat_clock_value = sat_clock.Galileo;
elseif use_GNSS == "QZSS"
    pseudorange_value = Data_OBS.pseudorange(:,121:129);
    sat_pos_value = sat_pos.QZSS;
    sat_clock_value = sat_clock.QZSS;
elseif use_GNSS == "All"
    pseudorange_value = Data_OBS.pseudorange;
    sat_pos_value = sat_pos.All;
    sat_clock_value = sat_clock.All;
%     sat_pos_value = Satpos;
end

user_pose_init = [0 0 0 0]';
error_threshold = 10e-3;
elev_mask = 15;
rms_error_ecef = zeros(epoch,3);
useful_sat_total = zeros(epoch,18);

for e = 1:epoch
    e
    resiErr = 10;

    user_pose = [0 0 0 0]';

    useful_sat_pse = find(isfinite(pseudorange_value(e,:)));
    useful_sat_ob = zeros(1,10);
    sat_use = find(isfinite(sat_pos_value(e,:)));

    j = 1;
    for i = 1:length(sat_use)
        if mod(sat_use(i),3) == 0
            useful_sat_ob(j) = sat_use(i)/3;
            j = j + 1;
        end
    end

    stop = 0;

    while(error_threshold < resiErr)
        useful_sat = intersect(useful_sat_pse,useful_sat_ob);
        useful_sat_num = length(useful_sat);

        x = zeros(useful_sat_num,3);
        rho = zeros(useful_sat_num,1);
        e_user = zeros(useful_sat_num,3);
        H = zeros(useful_sat_num,4);
        y = zeros(useful_sat_num,1);
       
        stop = stop + 1;
        if stop > 50
            break;
        end

        i = 1;
        while (i<=useful_sat_num)
            if useful_sat(i) >= 1 && useful_sat(i) <= 37
                sat_sys = "GPS";
            elseif useful_sat(i) >= 38 && useful_sat(i) <= 61
                sat_sys = "GLONASS";
            elseif useful_sat(i) >= 62 && useful_sat(i) <= 90
                sat_sys = "BeiDou";
            elseif useful_sat(i) >= 91 && useful_sat(i) <= 120
                sat_sys = "Galileo";
            else
                sat_sys = "QZSS";
            end

            x = sat_pos_value(e,3*useful_sat(i)-2:3*useful_sat(i));
            rho = pseudorange_value(e,useful_sat(i));
            b = sat_clock_value(e,useful_sat(i)); % sat clock bias

            % elevation angle cutoff 
            [elev,~] = calelevation(x,user_pose(1:3)');

            if stop > 1 && abs(elev) < 15
                useful_sat(i) = [];                
                useful_sat_num = useful_sat_num - 1;
                continue;
            end

            % Pseudorange correction
            corr_rho =   rho + speed_of_light*b;
            
            % Satellite Flight Correction
            dtflight = (corr_rho - user_pose(4))/speed_of_light + b;
            corr_satpos = RotatCorr(x,dtflight,sat_sys);

            % LS Componant
            d = sqrt((corr_satpos(1) - user_pose(1))^2 + (corr_satpos(2) - user_pose(2))^2 + (corr_satpos(3) - user_pose(3))^2);
            e_user(i,:) = (corr_satpos - user_pose(1:3)')/d;
            H(i,:) = [e_user(i,:) -1];
            y(i,1) = dot(corr_satpos,e_user(i,:)) - corr_rho;

            
            % test_sat_pos(i,:) = corr_satpos;
            % test_pseudo(i,1) = corr_rho;
            i = i + 1;

        end

        user_pose_old = user_pose;
        user_pose = (H'*H)\H'*y;

        resiErr = norm(user_pose_old - user_pose);
    end
    % recPos(e,:) = receiverposition(test_pseudo,test_sat_pos);
    cal_user_pose(e,:) = user_pose;
    [cal_user_pose_geodetic(e,1), cal_user_pose_geodetic(e,2), cal_user_pose_geodetic(e,3)] = ecef2geodetic(wgs84,user_pose(1),user_pose(2),user_pose(3));
    cal_user_pose_lla(e,:) = ecef2lla(user_pose(1:3)');

    if useful_sat_num ~= 0
        useful_sat_total(e,1:useful_sat_num) = useful_sat;
    end

    rms_error_ecef(e,:) = rmse(cal_user_pose(:,1:3), ref_pos);
end

%% error
% rms_error_ecef = rmse(cal_user_pose(:,1:3), ref_pos);

for i = 1:epoch
    [error_lla(i,1), ~]= lldistm(cal_user_pose_lla(i,1:2),[ref_pos_geodetic(1) ref_pos_geodetic(2)]);
end

rms_error_lla_hori = rms(error_lla);
rms_error_lla_height = rms(cal_user_pose_lla(:,3) - ref_pos_geodetic(3));


% %%
% figure
% % plot3(sat_pos(:,3*PRN(2)-2),sat_pos(:,3*PRN(2)-1),sat_pos(:,3*PRN(2)))
% for i = 1:useful_sat_num
%     plot3(sat_pos.GPS(:,3*useful_sat(i)-2)/1000,sat_pos.GPS(:,3*useful_sat(i)-1)/1000,sat_pos.GPS(:,3*useful_sat(i))/1000)
%     hold on
% end
% run('Earth_plot_nf');
% axis equal
% grid on


%%
figure
geoscatter(cal_user_pose_lla(:,1),cal_user_pose_lla(:,2),"filled"); geobasemap satellite
hold on
geoscatter(ref_pos_geodetic(1),ref_pos_geodetic(2), 'o', 'LineWidth', 5)
legend('Positioning', 'True')
% figure
% geoscatter(cal_user_pose_geodetic(:,1),cal_user_pose_geodetic(:,2),"filled")

% figure
% plot3(cal_user_pose(:,1),cal_user_pose(:,2),cal_user_pose(:,3), '.')
% hold on
% plot3(ref_pos(1),ref_pos(2),ref_pos(3), 'o')
% hold on
% plot3(0,0,0,'o','LineWidth',3)
% grid on

% uif = uifigure;
% g = geoglobe(uif);
% geoplot3(g,cal_user_pose_geodetic(:,1),cal_user_pose_geodetic(:,2),cal_user_pose_geodetic(:,3))
% 

figure
tiledlayout(2,1)
nexttile
plot(rms_error_ecef(:,1),'LineWidth',3)
hold on
plot(rms_error_ecef(:,2),'LineWidth',3)
hold on
plot(rms_error_ecef(:,3),'LineWidth',3)
% ylim([0 30])
xlim([0 epoch])
grid on
xlabel('epoch')
ylabel('rms error(m)')
legend('x', 'y', 'z')
title('ECEF RMS Error(m)')

nexttile
plot(error_lla,'LineWidth',3)
hold on
plot(ref_pos_geodetic(3) - cal_user_pose_lla(:,3),'LineWidth',3)
grid on
xlabel('epoch')
xlim([0 epoch])
ylabel('range error(m)')
legend('horizontal', 'height')
title('LLA Range Error(m)')
