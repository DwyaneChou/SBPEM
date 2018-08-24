% Predict Correct
function [tau_n,U,V,Z,LU0,LV0,LZ0] = PC(pass,MESH,STATE,time_step)

U0 = STATE.U;
V0 = STATE.V;
Z0 = STATE.Z;

% F0
[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);

LU0  = LU;
LV0  = LV;
LZ0  = LZ;
LU = -LU;
LV = -LV;
LZ = -LZ;

% F1
STATE.U = U0 + 0.5*time_step*LU;
STATE.V = V0 + 0.5*time_step*LV;
STATE.Z = Z0 + 0.5*time_step*LZ;

[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
                                  
LU1 = -LU;
LV1 = -LV;
LZ1 = -LZ;

% F2
STATE.U = U0 + 0.5*time_step*LU1;
STATE.V = V0 + 0.5*time_step*LV1;
STATE.Z = Z0 + 0.5*time_step*LZ1;

[LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
                                  
LU2 = -LU;
LV2 = -LV;
LZ2 = -LZ;

LF1LF2    = inner_product(MESH,LU1,LV1,LZ1,LU2,LV2,LZ2);
LF2_norm2 = inner_product(MESH,LU2,LV2,LZ2,LU2,LV2,LZ2);

beta_n = LF1LF2/LF2_norm2;

if LF2_norm2 == 0
    tau_n = 0;
else
    tau_n = beta_n*time_step;
end

U = U0 + tau_n*LU2;
V = V0 + tau_n*LV2;
Z = Z0 + tau_n*LZ2;