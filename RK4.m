% 4th order Runge-Kutta
function [tau_n,U,V,Z,LU0,LV0,LZ0] = RK4(pass,MESH,STATE,time_step)

U0 = STATE.U;
V0 = STATE.V;
Z0 = STATE.Z;

% K1
[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
LU0  = LU;
LV0  = LV;
LZ0  = LZ;
U_k1 = -LU;
V_k1 = -LV;
Z_k1 = -LZ;

% K2
STATE.U = U0 + 0.5*time_step*U_k1;
STATE.V = V0 + 0.5*time_step*V_k1;
STATE.Z = Z0 + 0.5*time_step*Z_k1;

[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
                                  
U_k2 = -LU;
V_k2 = -LV;
Z_k2 = -LZ;

% K3
STATE.U = U0 + 0.5*time_step*U_k2;
STATE.V = V0 + 0.5*time_step*V_k2;
STATE.Z = Z0 + 0.5*time_step*Z_k2;

[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
                                  
U_k3 = -LU;
V_k3 = -LV;
Z_k3 = -LZ;

% K4
STATE.U = U0 + time_step*U_k3;
STATE.V = V0 + time_step*V_k3;
STATE.Z = Z0 + time_step*Z_k3;

[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
                                  
U_k4  = -LU;
V_k4  = -LV;
Z_k4  = -LZ;

phi4_U = (U_k1 + 2.0*U_k2 + 2.0*U_k3 + U_k4)/6.0;
phi4_V = (V_k1 + 2.0*V_k2 + 2.0*V_k3 + V_k4)/6.0;
phi4_Z = (Z_k1 + 2.0*Z_k2 + 2.0*Z_k3 + Z_k4)/6.0;

phi4_norm2 = inner_product(MESH,phi4_U,phi4_V,phi4_Z,phi4_U,phi4_V,phi4_Z);
        
R1R2       = inner_product(MESH,U_k1,V_k1,Z_k1,U_k2,V_k2,Z_k2);
        
R2R3       = inner_product(MESH,U_k2,V_k2,Z_k2,U_k3,V_k3,Z_k3);
        
R3R4       = inner_product(MESH,U_k3,V_k3,Z_k3,U_k4,V_k4,Z_k4);

beta_n = 1.0/(3.0*phi4_norm2).*(R1R2+R2R3+R3R4);
if phi4_norm2 == 0
    tau_n  = 0;
else
    tau_n  = beta_n*time_step;
end

U  = U0+tau_n*phi4_U;
V  = V0+tau_n*phi4_V;
Z  = Z0+tau_n*phi4_Z;
