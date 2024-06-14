Re = 6378.137;
load('topo.mat','topo','topomap1');
phi=linspace(0,pi,30);
theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);
xsphere=Re*sin(phi).*cos(theta);
ysphere=Re*sin(phi).*sin(theta);
zsphere=Re*cos(phi);
mhndl1=mesh(xsphere,ysphere,zsphere);
h = surface(xsphere,ysphere,zsphere,'FaceColor','texture','CData',topo);
axis equal
hold on;