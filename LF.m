% Leap Frog
function [tau_n,U_np1,V_np1,Z_np1,LU,LV,LZ] = LF(pass,time_step,nt,...
                                                 U,V,Z,U_nm1,V_nm1,Z_nm1,LU_nm1,LV_nm1,LZ_nm1,...
                                                 dlambda,dtheta,a,Omega,g,...
                                                 lat_u,lat_v,lat_z,...
                                                 nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                 coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y)

if nt ==1
    [tau_n,U_np1,V_np1,Z_np1,LU,LV,LZ] = RK4(pass,time_step,U,V,Z,dlambda,dtheta,a,Omega,g,...
                                             lat_u,lat_v,lat_z,...
                                             nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                             coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
else
                                             
U0 = U;
V0 = V;
Z0 = Z;

cosU = cos(lat_u);
cosV = cos(lat_v);
cosZ = cos(lat_z);

% LF_n
[LU,LV,LZ] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);

% LF_np1
U = U_nm1 - 2.0*time_step*LU;
V = V_nm1 - 2.0*time_step*LV;
Z = Z_nm1 - 2.0*time_step*LZ;

[LU_np1,LV_np1,LZ_np1] = spatial_discrete(pass,U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                                          nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                          coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);

% R1, R2, psi
U_R1   = (U - U0)/time_step;
U_R2   = (U0 - U_nm1)/time_step;

V_R1   = (V - V0)/time_step;
V_R2   = (V0 - V_nm1)/time_step;

Z_R1   = (Z - Z0)/time_step;
Z_R2   = (Z0 - Z_nm1)/time_step;

U_psi  = -(5.0*LU_np1 + 8.0*LU - LU_nm1)/12.0;
V_psi  = -(5.0*LV_np1 + 8.0*LV - LV_nm1)/12.0;
Z_psi  = -(5.0*LZ_np1 + 8.0*LZ - LZ_nm1)/12.0;

LFR1   = sum(sum(LU_np1.*U_R1.*cosU))+sum(sum(LV_np1.*V_R1.*cosV))+sum(sum(LZ_np1.*Z_R1.*cosZ));
LFR2   = sum(sum(LU_nm1.*U_R2.*cosU))+sum(sum(LV_nm1.*V_R2.*cosV))+sum(sum(LZ_nm1.*Z_R2.*cosZ));
psi2   = sum(sum(U_psi.*U_psi.*cosU))+sum(sum(V_psi.*V_psi.*cosV))+sum(sum(Z_psi.*Z_psi.*cosZ));

tau_n  = -time_step/(6.0*psi2)*(5.0*LFR1+LFR2);

U_np1  = U0 + tau_n*U_psi;
V_np1  = V0 + tau_n*V_psi;
Z_np1  = Z0 + tau_n*Z_psi;

end