function [LU,LV,LZ] = slow_pass(U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                                 nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                 coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)
% lat/lon coef.
cosLatV           = cos(lat_v);

% Prepare stagger field
h                 = sqrt(Z);
hip1(1:nx_z-1,:)  = h(2:nx_z,:);
hip1(nx_z    ,:)  = h(1,:);
hjp1(:,1:ny_z-1)  = h(:,2:ny_z);
hjp1(:,ny_z    )  = 0;% For pole

hOnU              = 0.5*(h+hip1); % h on u grid
hOnV_temp         = 0.5*(h+hjp1); % h on v grid
hOnV              = hOnV_temp(:,1:ny_v);

u                 = U./hOnU;
v                 = V./hOnV;

uOnV                  = u(:,1:ny_v);
uOnV_im1(2:nx_v,:)    = uOnV(1:nx_v-1,:);
uOnV_im1(1     ,:)    = uOnV(nx_v,:);
uOnV_jp1(:,1:ny_v-1)  = uOnV(:,2:ny_v);
uOnV_jp1(:,ny_v   )   = 0;
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

Vip1(1:nx_v-1,:)    = V(2:nx_v,:);
Vip1(nx_v    ,:)    = V(1,:);
Vim1(2:nx_v  ,:)    = V(1:nx_v-1,:);
Vim1(1       ,:)    = V(nx_v    ,:);
Vjp1(:,1:ny_v-1)    = V(:,2:ny_v);
Vjp1(:,ny_v    )    = 0;% For pole
Vjm1(:,2:ny_v  )    = V(:,1:ny_v-1);
Vjm1(:,1       )    = 0;% For pole

uU                  = u.*U;
uU_ip1(1:nx_u-1,:)  = uU(2:nx_u,:);
uU_ip1(nx_u    ,:)  = uU(1,:);
uU_im1(2:nx_u  ,:)  = uU(1:nx_u-1,:);
uU_im1(1       ,:)  = uU(nx_u,:);

vcos                     = v.*cosLatV;
vcosOnU                  = vcos;
vcosOnU(:,ny_u)          = 0;
vcosOnU_ip1(1:nx_u-1,:)  = vcosOnU(2:nx_u,:);
vcosOnU_ip1(nx_u    ,:)  = vcosOnU_ip1(1,:);
vcosOnU_jm1(:,2:ny_u  )  = vcosOnU(:,1:ny_u-1);
vcosOnU_jm1(:,1       )  = 0;% For pole
vcosOnU_ip1jm1(:,2:ny_u) = vcosOnU_ip1(:,1:ny_u-1);
vcosOnU_ip1jm1(:,1     ) = 0;% For pole

vVcos                 = vcos.*V;
vVcos_jp1(:,1:ny_v-1) = vVcos(:,2:ny_v);
vVcos_jp1(:,ny_v    ) = 0;% For pole
vVcos_jm1(:,2:ny_v  ) = vVcos(:,1:ny_v-1);
vVcos_jm1(:,1       ) = 0;% For pole

% Advection
ADV_U_x  = u      .*(Uip1-Uim1) + uU_ip1    - uU_im1;

ADV_U_y  = (vcosOnU_ip1 + vcosOnU).*Ujp1 - (vcosOnU_ip1jm1 + vcosOnU_jm1).*Ujm1;

ADV_V_x  = (uOnV_jp1    + uOnV   ).*Vip1 - (uOnV_im1jp1    + uOnV_im1   ).*Vim1;

ADV_V_y  = vcos   .*(Vjp1-Vjm1) + vVcos_jp1 - vVcos_jm1;

% Construct L operator
LU1         = ADV_U_x;
LU2         = ADV_U_y;
LU2(:,1  )  = 0;
LU2(:,end)  = 0;
LU          = coefU_x.*LU1+coefU_y.*LU2;
LU(:,1  )   = 0;% For southern pole
LU(:,end)   = 0;% For northern pole

LV1         = ADV_V_x;
LV2         = ADV_V_y;
LV          = coefV_x.*LV1+coefV_y.*LV2;

LZ          = 0;