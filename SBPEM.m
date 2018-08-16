% Spherical Barotropic Primitive Eqation Model
% Written by Zhou Lilong Jun 22,2018

clc
clear

time_start = clock;

% Choose test case from 'RH' or 'IM'
test_case = 'IM';

% Choose the integral scheme
% 1 for CDO(Consisitent Dissipation Operator)
% 2 for 4th order Runge-Kutta
% 3 for Predict Correct
% 4 for leap frog
IntSch = 3;

% Define the grid resolution
MESH.dx = 2.0; % Degree
MESH.dy = 2.0; % Degree

% Define time(seconds)
time_step     = 30.0;
run_time      = 30*24*3600;

% Choose the split parameter
% split_scheme = 1 for CSP2
% split_scheme = 0 for no split
Split.split_scheme = 0;
Split.split_num    = 8;

% Define output
history_interval = 3600;
output_precision = 'NC_FLOAT';

if Split.split_scheme == 1
    Split.fast_dt      = time_step/Split.split_num;
    Split.slow_dt      = 0.5*time_step;
end

% Generate C-grid on sphere
MESH = genMesh(MESH);

max_int_step = MESH.a*cosd(90-MESH.dy)*MESH.dlambda/300;

if strcmp(test_case,'RH')
    % Initial fields with Rossby-Haurwitz Wave
    [u,v,STATE.Z,MESH.ghs] = Haurwitz(MESH);
elseif strcmp(test_case,'IM')
    % Initial fields with Isolated Mountain
    [u,v,STATE.Z,MESH.ghs] = Isolated_Mountain(MESH);
end

% IAP transformation
STATE = IAP(MESH,STATE,u,v);

total_energy0 = inner_product(MESH,STATE.U,STATE.V,(STATE.Z+MESH.ghs),STATE.U,STATE.V,(STATE.Z+MESH.ghs));
total_mass0   = sum(sum(STATE.Z.*MESH.cosLatZ));

% Compute the coefficient for L operator
MESH.coefU_x = 0.25./(MESH.a * MESH.cosLatU * MESH.dlambda); %coefU_x = 1/(2*a*cos(theta_u)*2dx)
MESH.coefU_y = 0.25./(MESH.a * MESH.cosLatU * MESH.dtheta ); %coefU_y = 1/(2*a*cos(theta_u)*2dy)
MESH.coefV_x = 0.25./(MESH.a * MESH.cosLatV * MESH.dlambda); %coefV_x = 1/(2*a*cos(theta_v)*2dx)
MESH.coefV_y = 0.25./(MESH.a * MESH.cosLatV * MESH.dtheta ); %coefV_x = 1/(2*a*cos(theta_v)*2dy)
MESH.coefZ_x = 0.5 ./(MESH.a * MESH.cosLatZ * MESH.dlambda); %coefZ_x = 1/(  a*cos(theta_z)*2dx)
MESH.coefZ_y = 0.5 ./(MESH.a * MESH.cosLatZ * MESH.dtheta ); %coefZ_x = 1/(  a*cos(theta_z)*2dy)

% Output the initial status
int_step_num = ceil(run_time/time_step);
output_count = 0;
output_num   = ceil(run_time/history_interval)+1;

output_netCDF(MESH,STATE,history_interval,output_count,output_precision)

STATE_old.U  = 0;
STATE_old.V  = 0;
STATE_old.Z  = 0;
TEND_old.LU  = 0;
TEND_old.LV  = 0;
TEND_old.LZ  = 0;
ti_start     = clock;
sum_it       = 0;
for it = 1:int_step_num
    % Set previous tend and status for leap-frog
    if it>1
        STATE_old.U = STATE.U;
        STATE_old.V = STATE.V;
        STATE_old.Z = STATE.Z;
        TEND_old.LU = TEND.LU;
        TEND_old.LV = TEND.LV;
        TEND_old.LZ = TEND.LZ;
        STATE.U     = STATE_new.U;
        STATE.V     = STATE_new.V;
        STATE.Z     = STATE_new.Z;
    end
    
    STATE.intStep = it;
    
    if Split.split_scheme>0
        [tau_n,STATE_new.U,STATE_new.V,STATE_new.Z,...
               TEND.LU,TEND.LV,TEND.LZ] = split_integrator(MESH,STATE,STATE_old,TEND_old,Split,IntSch);
    else
        pass = 0;
        [tau_n,STATE_new.U,STATE_new.V,STATE_new.Z,...
               TEND.LU,TEND.LV,TEND.LZ] = integrator(pass,MESH,STATE,STATE_old,TEND_old,IntSch,time_step);
    end
    
    integral_time = it*time_step;
    
    % Output
	if rem(integral_time,history_interval)==0 && integral_time>=history_interval
        ti_end = clock;
                
        output_count = output_count+1;
        output_netCDF(MESH,STATE,history_interval,output_count,output_precision)
        
        sum_it  = sum_it + etime(ti_end,ti_start);
        ave_it  = sum_it/(output_count+1);
        rest_it = ave_it*(output_num-1-output_count)/60;
                  
        total_energy       = inner_product(MESH,STATE.U,STATE.V,(STATE.Z+MESH.ghs),STATE.U,STATE.V,(STATE.Z+MESH.ghs));
        total_mass         = sum(sum(STATE.Z.*MESH.cosLatZ));
        energy_change_rate = (total_energy-total_energy0)/total_energy0; %ECR
        mass_change_rate   = (total_mass-total_mass0)/total_mass0;       %MCR
        
        disp(['Output        = ',num2str(output_count),'/',num2str(output_num-1)]);
        disp(['tau_n         = ',num2str(tau_n)])
        disp(['Total Energy  = ',num2str(total_energy)]);
        disp(['Total Mass    = ',num2str(total_mass)]);
        disp(['ECR           = ',num2str(energy_change_rate)]);
        disp(['MCR           = ',num2str(mass_change_rate)]);
        disp(['Integral time = ',num2str(etime(ti_end,ti_start))]);
        disp(['Rest Int time = ',num2str(rest_it)]);
        disp('                                     ');
        
        ti_start = clock;
	end
end

time_end = clock;
disp(['It took ',num2str(etime(time_end,time_start)),' seconds to run SBPEM'])