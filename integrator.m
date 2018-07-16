% Choose the integral scheme
% 1 for HDO(harmonious diffusion operator)
% 2 for 4th order Runge-Kutta
% 3 for Predict Correct
% 4 for leap frog
function [tau_n,U,V,Z,LU,LV,LZ] = integrator(pass,MESH,STATE,STATE_old,TEND_old,IntSch,time_step)
if IntSch==1
    [tau_n,U,V,Z,LU,LV,LZ] = HDO(pass,MESH,STATE,time_step);
elseif IntSch==2
    [tau_n,U,V,Z,LU,LV,LZ] = RK4(pass,MESH,STATE,time_step);
elseif IntSch==3
    [tau_n,U,V,Z,LU,LV,LZ] = PC (pass,MESH,STATE,time_step);
elseif IntSch==4
    [tau_n,U,V,Z,LU,LV,LZ] = LF (pass,MESH,STATE,STATE_old,TEND_old,time_step);
end