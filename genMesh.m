function MESH = genMesh(MESH)
dx = MESH.dx;
dy = MESH.dy;

% Define Constants
MESH.Omega = 7.292*10^-5;        %2.0 * pi / 86400.0; 
MESH.a     = 6371220.0;          %6371229.0
MESH.g     = 9.80616;


% Generate C-grid on sphere, z represents the geopotential height
longitude_z = 0  :dx:360-dx;
latitude_z  = -90:dy:90;

longitude_u = longitude_z + 0.5*dx;
latitude_u  = latitude_z;

longitude_v = longitude_z;
latitude_v  = latitude_z(:,1:end-1) + 0.5*dy;

[lat_u,lon_u] = meshgrid(latitude_u,longitude_u);
[lat_v,lon_v] = meshgrid(latitude_v,longitude_v);
[lat_z,lon_z] = meshgrid(latitude_z,longitude_z);

% Get grid size
MESH.nx_u = size(lon_u,1);
MESH.ny_u = size(lon_u,2);
MESH.nx_v = size(lon_v,1);
MESH.ny_v = size(lon_v,2);
MESH.nx_z = size(lon_z,1);
MESH.ny_z = size(lon_z,2);

% Convert longitude/latitude from degree to radian
d2r        = pi/180.0;
MESH.lon_u = lon_u*d2r;
MESH.lat_u = lat_u*d2r;
MESH.lon_v = lon_v*d2r;
MESH.lat_v = lat_v*d2r;
MESH.lon_z = lon_z*d2r;
MESH.lat_z = lat_z*d2r;

MESH.dlambda = dx*d2r;
MESH.dtheta  = dy*d2r;

% lat/lon coef.
MESH.sinLatU                       = sin(MESH.lat_u);
MESH.sinLatU(:,1  )                = -1;
MESH.sinLatU(:,end)                = 1;
MESH.sinLatU_jp1(:,1:MESH.ny_u-1)  = MESH.sinLatU(:,2:MESH.ny_u);
MESH.sinLatU_jp1(:,MESH.ny_u    )  = 1;% For pole

MESH.sinLatV                       = sin(MESH.lat_v);

MESH.sinLatZ                       = sin(MESH.lat_z);
MESH.sinLatZ(:,1  )                = -1;
MESH.sinLatZ(:,end)                = 1;

MESH.tanLatU                       = tan(MESH.lat_u);
MESH.tanLatU(:,1  )                = 0;% For pole
MESH.tanLatU(:,end)                = 0;% For pole
MESH.tanLatU_jp1(:,1:MESH.ny_u-1)  = MESH.tanLatU(:,2:MESH.ny_u);
MESH.tanLatU_jp1(:,MESH.ny_u    )  = 0;% For pole

MESH.cosLatU                       = cos(MESH.lat_u);
MESH.cosLatV                       = cos(MESH.lat_v);
MESH.cosLatZ                       = cos(MESH.lat_z);

MESH.cosLatVOnU                    = MESH.cosLatV;
MESH.cosLatVOnU(:,MESH.ny_u)       = 0;% For pole
MESH.cosLatVOnU_jm1(:,2:MESH.ny_u) = MESH.cosLatVOnU(:,1:MESH.ny_u-1);
MESH.cosLatVOnU_jm1(:,1     )      = 0;% For pole

MESH.cosLatZOnV                    = MESH.cosLatZ(:,1:MESH.ny_v);
MESH.cosLatZOnV_jp1(:,1:MESH.ny_v) = MESH.cosLatZ(:,2:MESH.ny_z);

% Reset pole
MESH.cosLatZ(:,1  )                = 0.25*MESH.cosLatV(:,1);
MESH.cosLatZ(:,end)                = 0.25*MESH.cosLatV(:,end);
MESH.cosLatU(:,1  )                = 0.25*MESH.cosLatV(:,1);
MESH.cosLatU(:,end)                = 0.25*MESH.cosLatV(:,end);

% Plot V grid
% [x,y,z] = sph2cart(MESH.lon_v,MESH.lat_v,MESH.a);
% surf(x,y,z);