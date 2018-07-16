function [LU,LV,LZ] = fast_pass(STATE,MESH)

U = STATE.U;
V = STATE.V;
Z = STATE.Z;

Omega  = MESH.Omega;
a      = MESH.a;
dtheta = MESH.dtheta;

nx_u = MESH.nx_u;
ny_u = MESH.ny_u;
nx_v = MESH.nx_v;
ny_v = MESH.ny_v;
nx_z = MESH.nx_z;
ny_z = MESH.ny_z;

% Prepare stagger field
h                 = sqrt(Z);
hip1(1:nx_z-1,:)  = h(2:nx_z,:);
hip1(nx_z    ,:)  = h(1,:);
him1(2:nx_z  ,:)  = h(1:nx_z-1,:);
him1(1       ,:)  = h(nx_z,:);
hjp1(:,1:ny_z-1)  = h(:,2:ny_z);
hjp1(:,ny_z    )  = 0;% For pole
hjm1(:,2:ny_z  )  = h(:,1:ny_z-1);
hjm1(:,1       )  = 0;% For pole

hOnU              = 0.5*(h+hip1); % h on u grid
hOnV_temp         = 0.5*(h+hjp1); % h on v grid
hOnV              = hOnV_temp(:,1:ny_v);

u                 = U./hOnU;

Zip1(1:nx_z-1,:)  = Z(2:nx_z,:);
Zip1(nx_z    ,:)  = Z(1,:);

ZOnV                 = Z(:,1:ny_v);
ZOnV_jp1(:,1:ny_v-1) = ZOnV(:,2:ny_v);
ZOnV_jp1(:,ny_v    ) = Z(:,ny_z);

Uim1(2:nx_u  ,:)     = U(1:nx_u-1,:);
Uim1(1       ,:)     = U(nx_u,:);
Ujp1(:,1:ny_u-1)     = U(:,2:ny_u);
Ujp1(:,ny_u    )     = 0;% For pole
Uim1jp1(2:nx_u,:)    = Ujp1(1:nx_u-1,:);
Uim1jp1(1     ,:)    = Ujp1(nx_u,:);

Vip1(1:nx_v-1,:)    = V(2:nx_v,:);
Vip1(nx_v    ,:)    = V(1,:);
Vjm1(:,2:ny_v  )    = V(:,1:ny_v-1);
Vjm1(:,1       )    = 0;% For pole
Vip1jm1(1:nx_v-1,:) = Vjm1(2:nx_v,:);
Vip1jm1(nx_v    ,:) = Vjm1(1,:);

uU                  = u.*U;
uU_im1(2:nx_u  ,:)  = uU(1:nx_u-1,:);
uU_im1(1       ,:)  = uU(nx_u,:);
uU_jp1(:,1:ny_u-1)  = uU(:,2:ny_u);
uU_jp1(:,ny_u    )  = 0;% For pole
uU_im1jp1(2:nx_u,:) = uU_jp1(1:nx_u-1,:);
uU_im1jp1(1     ,:) = uU_jp1(nx_u,:);

VcosOnZ_temp          = V.*MESH.cosLatV;
VcosOnZ               = VcosOnZ_temp;
VcosOnZ(:,ny_z)       = 0;% For pole
VcosOnZ_jm1(:,2:ny_z) = VcosOnZ(:,1:ny_z-1);
VcosOnZ_jm1(:,1     ) = 0;% For pole

% Flux
FLUX_Z_x = (hip1+h).*U       - (h+him1).*Uim1;
FLUX_Z_y = (hjp1+h).*VcosOnZ - (h+hjm1).*VcosOnZ_jm1;

% Pressure Gradient Force
PGF_U             = 4.0*MESH.coefU_x.*hOnU.*(Zip1-Z);
PGF_U(:,1  )      = 0;
PGF_U(:,end)      = 0;

PGF_V_temp        = hOnV.*(ZOnV_jp1-ZOnV);
PGF_V             = PGF_V_temp/(a*dtheta);

% Colioris Force
C1                      = MESH.cosLatVOnU    ./MESH.cosLatZ;
C2                      = MESH.cosLatVOnU_jm1./MESH.cosLatZ;
VOnU_part1_temp         = Vip1+V;
VOnU_part1_temp(:,ny_u) = 0;
VOnU_part1              = VOnU_part1_temp.*C1;
VOnU_part2_temp         = Vip1jm1+Vjm1;
VOnU_part2_temp(:,ny_u) = 0;
VOnU_part2              = VOnU_part2_temp.*C2;
VOnU                    = 0.25*(VOnU_part1+VOnU_part2);

fv                 = 2.0*Omega*MESH.sinLatU.*VOnU;

fu_temp            = 0.5*Omega*((Ujp1 + Uim1jp1).*MESH.sinLatU_jp1...
                                    +(U    + Uim1   ).*MESH.sinLatU    ); % 2*Omega*sin(theta)U/4
fu                 = fu_temp(:,1:ny_v);

% Curvature Term
CV                 = u/a.*MESH.tanLatU.*VOnU; % u/a*tan(theta)*V_On_U

CU_part1           = (uU_jp1 + uU_im1jp1).*MESH.tanLatU_jp1;
CU_part2           = (uU     + uU_im1   ).*MESH.tanLatU;
CU                 = 0.25*(CU_part1(:,1:ny_v)+CU_part2(:,1:ny_v))/a;

% Construct L operator
LU          = PGF_U-fv-CV;
LU(:,1  )   = 0;% For southern pole
LU(:,end)   = 0;% For northern pole

LV          = PGF_V+fu+CU;

LZ1         = MESH.coefZ_x.*FLUX_Z_x;
LZ2         = MESH.coefZ_y.*FLUX_Z_y;
LZ2(:,1  )  =  sum((hjp1(:,1  ) + h   (:,1  )).*V(:,1  ).*MESH.cosLatV(:,1  ))/nx_z.*MESH.coefZ_y(:,1  );% Southern Pole
LZ2(:,end)  = -sum((h   (:,end) + hjm1(:,end)).*V(:,end).*MESH.cosLatV(:,end))/nx_z.*MESH.coefZ_y(:,end);% Northern Pole
LZ          = LZ1+LZ2;
