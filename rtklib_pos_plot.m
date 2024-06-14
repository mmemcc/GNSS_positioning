clear all
ref_pos = 1.0e+06 * [-3.052949583000000   4.042646515000000   3.862274494000000];
ref_pos_geodetic = ecef2lla(ref_pos);
pos = readmatrix("Data\pos.csv");

pos(:,1) = [];
pos(1:14,:) = [];

figure
geoscatter(pos(:,1),pos(:,2),"filled");  geobasemap satellite
hold on
geoscatter(ref_pos_geodetic(1),ref_pos_geodetic(2), 'o', 'LineWidth', 5)
legend('Positioning', 'True')