function [LU1,LV1,LZ1,BU,BV,BZ] = B_operator(pass,STATE,MESH,time_step)

U0 = STATE.U;
V0 = STATE.V;
Z0 = STATE.Z;

% K1
[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
LU1  = LU;
LV1  = LV;
LZ1  = LZ;

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

BU     = 2.0 / time_step * (phi4_U - U_k1);
BV     = 2.0 / time_step * (phi4_V - V_k1);
BZ     = 2.0 / time_step * (phi4_Z - Z_k1);