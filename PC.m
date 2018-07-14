% Predict Correct
function [tau_n,U_np1,V_np1,Z_np1,LU0,LV0,LZ0] = PC(pass,time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                                    lat_u,lat_v,lat_z,...
                                                    nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                    coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)

U0 = U;
V0 = V;
Z0 = Z;

cosU = cos(lat_u);
cosV = cos(lat_v);
cosZ = cos(lat_z);

% F0
[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
LU0  = LU;
LV0  = LV;
LZ0  = LZ;
LU = -LU;
LV = -LV;
LZ = -LZ;

% F1
U  = U0 + 0.5*time_step*LU;
V  = V0 + 0.5*time_step*LV;
Z  = Z0 + 0.5*time_step*LZ;

[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
LU1 = -LU;
LV1 = -LV;
LZ1 = -LZ;

% F2
U   = U0 + 0.5*time_step*LU1;
V   = V0 + 0.5*time_step*LV1;
Z   = Z0 + 0.5*time_step*LZ1;

[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
LU2 = -LU;
LV2 = -LV;
LZ2 = -LZ;

LF1LF2   = sum(sum(LU1.*LU2.*cosU)) + sum(sum(LV1.*LV2.*cosV)) +sum(sum(LZ1.*LZ2.*cosZ));
LF2_nor2 = sum(sum(LU2.*LU2.*cosU)) + sum(sum(LV2.*LV2.*cosV)) +sum(sum(LZ2.*LZ2.*cosZ));

beta_n = LF1LF2/LF2_nor2;
tau_n  = beta_n*time_step;

U_np1  = U0 + tau_n*LU2;
V_np1  = V0 + tau_n*LV2;
Z_np1  = Z0 + tau_n*LZ2;