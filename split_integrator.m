function [tau_n,U,V,Z,LU,LV,LZ] = split_integrator(MESH,STATE,STATE_old,TEND_old,Split,IntSch)

fast_pass = 1;
slow_pass = 2;

if Split.split_scheme == 1
    
    dt = Split.slow_dt;
    [~,STATE.U,STATE.V,STATE.Z,...
       TEND.LU,TEND.LV,TEND.LZ] = integrator(slow_pass,MESH,STATE,STATE_old,TEND_old,IntSch,dt);
    
    for i = 1:Split.split_num
        if IntSch == 4
            STATE_old.U  = STATE.U;
            STATE_old.V  = STATE.V;
            STATE_old.Z  = STATE.Z;
            TEND_old.LU  = TEND.LU;
            TEND_old.LV  = TEND.LV;
            TEND_old.LZ  = TEND.LZ;
        end
        
        dt = Split.fast_dt;
        [~,STATE.U,STATE.V,STATE.Z,...
           TEND.LU,TEND.LV,TEND.LZ] = integrator(fast_pass,MESH,STATE,STATE_old,TEND_old,IntSch,dt);
    end
    
    if IntSch == 4
        STATE_old.U  = STATE.U;
        STATE_old.V  = STATE.V;
        STATE_old.Z  = STATE.Z;
        TEND_old.LU  = TEND.LU;
        TEND_old.LV  = TEND.LV;
        TEND_old.LZ  = TEND.LZ;
    end
    
    dt = Split.slow_dt;
    [tau_n,U,V,Z,LU,LV,LZ] = integrator(slow_pass,MESH,STATE,STATE_old,TEND_old,IntSch,dt);
end

tau_n = tau_n*2;