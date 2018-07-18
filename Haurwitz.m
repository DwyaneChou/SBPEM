function [u,v,Z] = Haurwitz(MESH)

MESH.cosLatZ(:,1  ) = 0;
MESH.cosLatZ(:,end) = 0;
MESH.cosLatU(:,1  ) = 0;
MESH.cosLatU(:,end) = 0;

a     = MESH.a;
Omega = MESH.Omega;
g     = MESH.g;
lon_u = MESH.lon_u;
lon_v = MESH.lon_v;
lon_z = MESH.lon_z;

% Set Constants
omega = 7.848*10^-6;
K     = omega;
R     = 4.0;
h0    = 8.0*10^3;

% Initial fields with Rossby-Haurwitz Wave, reference:
% "A Standard Test Set for Numerical Approximations to the Shallow
% Water Equations in Spherical Geometry"
u1  = MESH.cosLatU;
u2  = R*MESH.cosLatU.^(R-1).*MESH.sinLatU.^2.*cos(R*lon_u);
u3  =   MESH.cosLatU.^(R+1).*cos(R*lon_u);
u   = a*omega*(u1+u2-u3);

v   = -a*K*R*MESH.cosLatV.^(R-1).*MESH.sinLatV.*sin(R*lon_v);

A1  = omega*0.5*(2*Omega+omega)*MESH.cosLatZ.^2;
Ac  = 0.25*K^2;
A21 = (R+1.0).*MESH.cosLatZ.^(2.0*R+2);
A22 = (2.0*R^2-R-2)*MESH.cosLatZ.^(2.0*R);
A23 = 2.0*R^2.*MESH.cosLatZ.^(2.0*R-2);
A   = A1+Ac.*(A21+A22-A23);

Bc  = 2.*(Omega+omega)*K/((R+1)*(R+2)).*MESH.cosLatZ.^R;
B1  = R^2+2.0*R+2.0;
B2  = (R+1.0)^2.*MESH.cosLatZ.^2;
B   = Bc.*(B1-B2);

Cc  = 0.25*K^2.*MESH.cosLatZ.^(2.0*R);
C1  = (R+1.0)*MESH.cosLatZ.^2;
C2  = R+2.0;
C   = Cc.*(C1-C2);

Z  = g*h0+a^2*(A + B.*cos(R*lon_z) + C.*cos(2.0*R*lon_z));