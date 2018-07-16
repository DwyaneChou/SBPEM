function [LU,LV,LZ] = L_operator(STATE,MESH)

Omega = MESH.Omega;
a     = MESH.a;

U = STATE.U;
V = STATE.V;
Z = STATE.Z;

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
v                 = V./hOnV;

Zip1(1:nx_z-1,:)  = Z(2:nx_z,:);
Zip1(nx_z    ,:)  = Z(1,:);

ZOnV                  = Z(:,1:ny_v);
ZOnV_jp1(:,1:ny_v-1)  = ZOnV(:,2:ny_v);
ZOnV_jp1(:,ny_v    )  = Z(:,ny_z);

uOnV                  = u(:,1:ny_v);
uOnV_im1(2:nx_v,:)    = uOnV(1:nx_v-1,:);
uOnV_im1(1     ,:)    = uOnV(nx_v,:);
uOnV_jp1(:,1:ny_v-1)  = uOnV(:,2:ny_v);
uOnV_jp1(:,ny_v    )  = u(:,ny_u);
uOnV_im1jp1(2:nx_v,:) = uOnV_jp1(1:nx_v-1,:);
uOnV_im1jp1(1     ,:) = uOnV_jp1(nx_v,:);

Uip1(1:nx_u-1,:)      = U(2:nx_u,:);
Uip1(nx_u    ,:)      = U(1,:);
Uim1(2:nx_u  ,:)      = U(1:nx_u-1,:);
Uim1(1       ,:)      = U(nx_u,:);
Ujp1(:,1:ny_u-1)      = U(:,2:ny_u);
Ujp1(:,ny_u    )      = 0;% For pole
Ujm1(:,2:ny_u  )      = U(:,1:ny_u-1);
Ujm1(:,1       )      = 0;% For pole
Uim1jp1(2:nx_u,:)     = Ujp1(1:nx_u-1,:);
Uim1jp1(1     ,:)     = Ujp1(nx_u,:);

Vip1(1:nx_v-1,:)      = V(2:nx_v,:);
Vip1(nx_v    ,:)      = V(1,:);
Vim1(2:nx_v  ,:)      = V(1:nx_v-1,:);
Vim1(1       ,:)      = V(nx_v    ,:);
Vjp1(:,1:ny_v-1)      = V(:,2:ny_v);
Vjp1(:,ny_v    )      = 0;% For pole
Vjm1(:,2:ny_v  )      = V(:,1:ny_v-1);
Vjm1(:,1       )      = 0;% For pole
Vip1jm1(1:nx_v-1,:)   = Vjm1(2:nx_v,:);
Vip1jm1(nx_v    ,:)   = Vjm1(1,:);

uU                    = u.*U;
uU_ip1(1:nx_u-1,:)    = uU(2:nx_u,:);
uU_ip1(nx_u    ,:)    = uU(1,:);
uU_im1(2:nx_u  ,:)    = uU(1:nx_u-1,:);
uU_im1(1       ,:)    = uU(nx_u,:);
uU_jp1(:,1:ny_u-1)    = uU(:,2:ny_u);
uU_jp1(:,ny_u    )    = 0;% For pole
uU_im1jp1(2:nx_u,:)   = uU_jp1(1:nx_u-1,:);
uU_im1jp1(1     ,:)   = uU_jp1(nx_u,:);

vcos                     = v.*MESH.cosLatV;
vcosOnU                  = vcos;
vcosOnU(:,ny_u)          = 0;
vcosOnU_ip1(1:nx_u-1,:)  = vcosOnU(2:nx_u,:);
vcosOnU_ip1(nx_u    ,:)  = vcosOnU_ip1(1,:);
vcosOnU_jm1(:,2:ny_u  )  = vcosOnU(:,1:ny_u-1);
vcosOnU_jm1(:,1       )  = 0;% For pole
vcosOnU_ip1jm1(:,2:ny_u) = vcosOnU_ip1(:,1:ny_u-1);
vcosOnU_ip1jm1(:,1     ) = 0;% For pole

vVcos                    = vcos.*V;
vVcos_jp1(:,1:ny_v-1)    = vVcos(:,2:ny_v);
vVcos_jp1(:,ny_v    )    = 0;% For pole
vVcos_jm1(:,2:ny_v  )    = vVcos(:,1:ny_v-1);
vVcos_jm1(:,1       )    = 0;% For pole

VcosOnZ                  = V.*MESH.cosLatV;
VcosOnZ(:,ny_z)          = 0;% For pole
VcosOnZ_jm1(:,2:ny_z)    = VcosOnZ(:,1:ny_z-1);
VcosOnZ_jm1(:,1     )    = 0;% For pole

% Advection
ADV_U_x  = u    .*(Uip1-Uim1) + uU_ip1    - uU_im1;

ADV_U_y  = (vcosOnU_ip1 + vcosOnU).*Ujp1 - (vcosOnU_ip1jm1 + vcosOnU_jm1).*Ujm1;

ADV_V_x  = (uOnV_jp1    + uOnV   ).*Vip1 - (uOnV_im1jp1    + uOnV_im1   ).*Vim1;

ADV_V_y  = vcos .*(Vjp1-Vjm1) + vVcos_jp1 - vVcos_jm1;

% Pressure Gradient Force
PGF_U             = 4.0*MESH.coefU_x.*hOnU.*(Zip1-Z);
PGF_U(:,1  )      = 0;
PGF_U(:,end)      = 0;

PGF_V_temp        = hOnV.*(ZOnV_jp1-ZOnV);
PGF_V             = PGF_V_temp/(MESH.a*MESH.dtheta);

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

% Flux
FLUX_Z_x           = (hip1+h).*U - (h+him1).*Uim1;

FLUX_Z_y_part1     =  (hjp1 + h   ).*VcosOnZ;
FLUX_Z_y_part2     = -(h    + hjm1).*VcosOnZ_jm1;
FLUX_Z_y           =  FLUX_Z_y_part1+FLUX_Z_y_part2;

% Construct L operator
LU1         = ADV_U_x;
LU2         = ADV_U_y;
LU2(:,1  )  = 0;
LU2(:,end)  = 0;
LU          = MESH.coefU_x.*LU1...
             +MESH.coefU_y.*LU2...
             +PGF_U-fv-CV;
LU(:,1  )   = 0;% For southern pole
LU(:,end)   = 0;% For northern pole

LV1         = ADV_V_x;
LV2         = ADV_V_y;
LV          = MESH.coefV_x.*LV1...
             +MESH.coefV_y.*LV2...
             +PGF_V+fu+CU;

LZ1         = MESH.coefZ_x.*FLUX_Z_x;
LZ2         = MESH.coefZ_y.*FLUX_Z_y;
LZ2(:,1  )  = mean(FLUX_Z_y_part1(:,1  ))*MESH.coefZ_y(:,1  );% Southern Pole
LZ2(:,end)  = mean(FLUX_Z_y_part2(:,end))*MESH.coefZ_y(:,end);% Northern Pole
LZ          = LZ1+LZ2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Symmetry Check %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% LUU         = LU.*U;
% LVV         = LV.*V;
% LZZ         = LZ.*Z;
% 
% figure
% pcolor(LUU)
% shading interp
% 
% figure
% pcolor(LVV)
% shading interp
% 
% figure
% pcolor(LZZ)
% shading interp
% 
% CFUU = fv.*U;
% CFVV = fu.*V;
% CTUU = CV.*U;
% CTVV = CU.*V;
% 
% LU1_U = sum(sum(LU1.*U));
% LU2_U = sum(sum(LU2.*U));
% LV1_V = sum(sum(LV1.*V));
% LV2_V = sum(sum(LV2.*V));
% PGUU  = sum(sum(PGF_U.*U.*MESH.cosLatU));
% LZ1Z  = sum(sum(LZ1.*Z.*MESH.cosLatZ));
% PGVV  = sum(sum(PGF_V.*V.*MESH.cosLatV));
% LZ2Z  = sum(sum(LZ2.*Z.*MESH.cosLatZ));
% fVU   = sum(sum(CFUU.*MESH.cosLatU));
% fUV   = sum(sum(CFVV.*MESH.cosLatV));
% CUV   = sum(sum(CTUU.*MESH.cosLatU));
% CVU   = sum(sum(CTVV.*MESH.cosLatV));
% LFF   = sum(sum(LUU.*MESH.cosLatU))+sum(sum(LVV.*MESH.cosLatV))+sum(sum(LZZ.*MESH.cosLatZ));
% 
% disp(['(LU1,U)     = ',num2str(LU1_U)])
% disp(['(LU2,U)     = ',num2str(LU2_U)])
% disp(['(LV1,V)     = ',num2str(LV1_V)])
% disp(['(LV2,V)     = ',num2str(LV2_V)])
% disp(['(PGUU+LZ1Z) = ',num2str(PGUU+LZ1Z)]);
% disp(['(PGVV+LZ2Z) = ',num2str(PGVV+LZ2Z)]);
% 
% disp(['(fVU+fUV)   = ',num2str(-fVU+fUV)]);
% disp(['(CUV+CVU)   = ',num2str(-CUV+CVU)]);
% 
% disp(['(LF,F)      = ',num2str(LFF)]);
% disp('                                  ');
% ;