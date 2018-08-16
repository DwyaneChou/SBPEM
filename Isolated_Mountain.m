function [u,v,Z,ghs] = Isolated_Mountain(MESH)

MESH.cosLatZ(:,1  ) = 0;
MESH.cosLatZ(:,end) = 0;
MESH.cosLatU(:,1  ) = 0;
MESH.cosLatU(:,end) = 0;

a     = MESH.a;
Omega = MESH.Omega;
g     = MESH.g;
lon_u = MESH.lon_u;
lat_u = MESH.lat_u;
lon_v = MESH.lon_v;
lon_z = MESH.lon_z;
lat_z = MESH.lat_z;

% Set Constants
hs0      = 2000;
R        = ones(size(lon_z))*pi/9;
labmda_c = 3*pi/2;
theta_c  = pi/6;
r        = sqrt(min(R.^2,(lon_z-labmda_c).^2+(lat_z-theta_c).^2));
ghs      = g*hs0*(1-r./R);

u0       = 20;%2*pi*a/(12*24*3600);
alpha    = 0; % 0 or 0.005 or pi/2-0.05 or pi/2
gh0      = 5960*g;%5400*g;

u        = u0*(cos(lat_u).*cos(alpha)+cos(lon_u).*sin(lat_u).*sin(alpha));
v        = -u0*sin(lon_v)*sin(alpha);
Z        = gh0 - (a*Omega*u0 + u0^2/2).*(-cos(lon_z).*cos(lat_z).*sin(alpha) + sin(lat_z).*cos(alpha)).^2;