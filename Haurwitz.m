function [u,v,Z] = Haurwitz(MESH)
a     = MESH.a;
Omega = MESH.Omega;
g     = MESH.g;
lon_u = MESH.lon_u;
lon_v = MESH.lon_v;
lon_z = MESH.lon_z;
lat_u = MESH.lat_u;
lat_v = MESH.lat_v;
lat_z = MESH.lat_z;

% Set Constants
omega = 7.848*10^-6;
K     = omega;
R     = 4.0;
h0    = 8.0*10^3;

% Initial fields with Rossby-Haurwitz Wave, reference:
% "A Standard Test Set for Numerical Approximations to the Shallow
% Water Equations in Spherical Geometry"
u1  = cos(lat_u);
u2  = R*cos(lat_u).^(R-1).*sin(lat_u).^2.*cos(R*lon_u);
u3  =   cos(lat_u).^(R+1).*cos(R*lon_u);
u   = a*omega*(u1+u2-u3);

v   = -a*K*R*cos(lat_v).^(R-1).*sin(lat_v).*sin(R*lon_v);

A1  = omega*0.5*(2*Omega+omega)*cos(lat_z).^2;
Ac  = 0.25*K^2*cos(lat_z).^(2.0*R);
A21 = (R+1.0).*cos(lat_z).^2;
A22 = 2.0*R^2-R-2;
A23 = 2.0*R^2.*cos(lat_z).^-2;
A   = A1+Ac.*(A21+A22-A23);

Bc  = 2.*(Omega+omega)*K/((R+1)*(R+2)).*cos(lat_z).^R;
B1  = R^2+2.0*R+2.0;
B2  = (R+1.0)^2.*cos(lat_z).^2;
B   = Bc.*(B1-B2);

Cc  = 0.25*K^2.*cos(lat_z).^(2.0*R);
C1  = (R+1.0)*cos(lat_z).^2;
C2  = R+2.0;
C   = Cc.*(C1-C2);

Z  = g*h0+a^2*(A + B.*cos(R*lon_z) + C.*cos(2.0*R*lon_z));