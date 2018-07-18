% Leap Frog
function [tau_n,U,V,Z,LU,LV,LZ] = LF(pass,MESH,STATE,STATE_old,TEND_old,time_step)

if STATE.intStep ==1
    [tau_n,U,V,Z,LU,LV,LZ] = RK4(pass,MESH,STATE,time_step);
else
    
    U0 = STATE.U;
    V0 = STATE.V;
    Z0 = STATE.Z;
    
    U_nm1  = STATE_old.U;
    V_nm1  = STATE_old.V;
    Z_nm1  = STATE_old.Z;
    
    LU_nm1 = TEND_old.LU;
    LV_nm1 = TEND_old.LV;
    LZ_nm1 = TEND_old.LZ;
    
    % LF_n
    [LU,LV,LZ] = spatial_discrete(pass,STATE,MESH);
    
    % LF_np1
    STATE.U = U_nm1 - 2.0*time_step*LU;
    STATE.V = V_nm1 - 2.0*time_step*LV;
    STATE.Z = Z_nm1 - 2.0*time_step*LZ;
    
    [LU_np1,LV_np1,LZ_np1] = spatial_discrete(pass,STATE,MESH);
    
    % R1, R2, psi
    U_R1   = (STATE.U - U0   )/time_step;
    U_R2   = (U0      - U_nm1)/time_step;
    
    V_R1   = (STATE.V - V0   )/time_step;
    V_R2   = (V0      - V_nm1)/time_step;
    
    Z_R1   = (STATE.Z - Z0   )/time_step;
    Z_R2   = (Z0      - Z_nm1)/time_step;
    
    U_psi  = -(5.0*LU_np1 + 8.0*LU - LU_nm1)/12.0;
    V_psi  = -(5.0*LV_np1 + 8.0*LV - LV_nm1)/12.0;
    Z_psi  = -(5.0*LZ_np1 + 8.0*LZ - LZ_nm1)/12.0;
    
    LFR1   = inner_product(MESH,LU_np1,LV_np1,LZ_np1,U_R1,V_R1,Z_R1);
    LFR2   = inner_product(MESH,LU_nm1,LV_nm1,LZ_nm1,U_R2,V_R2,Z_R2);
    psi2   = inner_product(MESH,U_psi,V_psi,Z_psi,U_psi,V_psi,Z_psi);
    
    tau_n  = -time_step/(6.0*psi2)*(5.0*LFR1+LFR2);
    
    U = U0 + tau_n*U_psi;
    V = V0 + tau_n*V_psi;
    Z = Z0 + tau_n*Z_psi;

end