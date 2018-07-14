% 4th order Runge-Kutta
function [tau_n,U_np1,V_np1,Z_np1,LU0,LV0,LZ0] = RK4(pass,time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                                     lat_u,lat_v,lat_z,...
                                                     nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                     coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)

U0 = U;
V0 = V;
Z0 = Z;

cosU = cos(lat_u);
cosV = cos(lat_v);
cosZ = cos(lat_z);

% K1
[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
LU0  = LU;
LV0  = LV;
LZ0  = LZ;
U_k1 = -LU;
V_k1 = -LV;
Z_k1 = -LZ;

% K2
U    = U0 + 0.5*time_step*U_k1;
V    = V0 + 0.5*time_step*V_k1;
Z    = Z0 + 0.5*time_step*Z_k1;

[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k2 = -LU;
V_k2 = -LV;
Z_k2 = -LZ;

% K3
U    = U0 + 0.5*time_step*U_k2;
V    = V0 + 0.5*time_step*V_k2;
Z    = Z0 + 0.5*time_step*Z_k2;

[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k3 = -LU;
V_k3 = -LV;
Z_k3 = -LZ;

% K4
U    = U0 + time_step*U_k3;
V    = V0 + time_step*V_k3;
Z    = Z0 + time_step*Z_k3;

[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
                                  
U_k4  = -LU;
V_k4  = -LV;
Z_k4  = -LZ;

phi4_U = (U_k1 + 2.0*U_k2 + 2.0*U_k3 + U_k4)/6.0;
phi4_V = (V_k1 + 2.0*V_k2 + 2.0*V_k3 + V_k4)/6.0;
phi4_Z = (Z_k1 + 2.0*Z_k2 + 2.0*Z_k3 + Z_k4)/6.0;

phi4_norm2 = sum(sum(phi4_U.*phi4_U.*cosU))+sum(sum(phi4_V.*phi4_V.*cosV))+sum(sum(phi4_Z.*phi4_Z.*cosZ));
R1R2       = sum(sum(U_k1.*U_k2.*cosU))+sum(sum(V_k1.*V_k2.*cosV))+sum(sum(Z_k1.*Z_k2.*cosZ));
R2R3       = sum(sum(U_k2.*U_k3.*cosU))+sum(sum(V_k2.*V_k3.*cosV))+sum(sum(Z_k2.*Z_k3.*cosZ));
R3R4       = sum(sum(U_k3.*U_k4.*cosU))+sum(sum(V_k3.*V_k4.*cosV))+sum(sum(Z_k3.*Z_k4.*cosZ));

beta_n = 1.0/(3.0*phi4_norm2).*(R1R2+R2R3+R3R4);
tau_n  = beta_n*time_step;

U_np1  = U0+tau_n*phi4_U;
V_np1  = V0+tau_n*phi4_V;
Z_np1  = Z0+tau_n*phi4_Z;
