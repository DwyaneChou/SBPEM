% Harmonious Diffusion Operator
function [tau_n,U_np1,V_np1,Z_np1,LU0,LV0,LZ0] = HDO(pass,MESH,STATE,time_step)

[LU,LV,LZ,BU,BV,BZ] = B_operator(pass,STATE,MESH,time_step);
LU0  = LU;
LV0  = LV;
LZ0  = LZ;

U = STATE.U;
V = STATE.V;
Z = STATE.Z;

% inner product
LF_norm2 = inner_product(MESH,LU,LV,LZ,LU,LV,LZ);
BFF      = inner_product(MESH,BU,BV,BZ,U ,V ,Z );
BFLF     = inner_product(MESH,BU,BV,BZ,LU,LV,LZ);
BF_norm2 = inner_product(MESH,BU,BV,BZ,BU,BV,BZ);

K1       = -LF_norm2/BFF;
K2       = BFLF/BFF;
K3       = -BF_norm2/BFF;

eps      = 0.5;
tau1     = 2.0 - K1 / eps;
tau2     = K2 + sqrt( K2^2 + eps * tau1 * K3 );
tau_n    = tau1/tau2;

U_np1    = U - LU * tau_n + eps * tau_n.^2 * BU;
V_np1    = V - LV * tau_n + eps * tau_n.^2 * BV;
Z_np1    = Z - LZ * tau_n + eps * tau_n.^2 * BZ;