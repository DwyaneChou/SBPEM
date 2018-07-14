% Spherical Barotropic Primitive Eqation Model
% Written by Zhou Lilong Jun 22,2018

clc
clear

time_start = clock;

% Choose the integral scheme
% 1 for HDO(harmonious diffusion operator)
% 2 for 4th order Runge-Kutta
% 3 for Predict Correct
% 4 for leap frog
IntSch = 2;

% Define Constants
Omega = 7.292*10^-5;
a     = 6371229.0;
g     = 9.8;

% Define the grid resolution
dx = 2.0; % Degree
dy = 2.0; % Degree

% Define time(seconds)
max_time_step = cosd(90-dy)*a/180.0*pi/1000*6;
time_step     = 18.0;
run_time      = 33*24*3600;

% Define output
history_interval = 3600;
output_precision = 'NC_FLOAT';

% Choose the split parameter
% split_scheme = 1 for CSP2
% split_scheme = 0 for no split
split_scheme = 0;
split_num    = 3;
fast_dt      = time_step/split_num;
slow_dt      = time_step;

% Generate C-grid on sphere, z represents the geopotential height
longitude_z = 0:dx:360-dx;
latitude_z  = -90:dx:90;

longitude_u = 0+0.5*dx:dx:360+0.5*dx-dx;
latitude_u  = -90:dx:90;

longitude_v = 0:dx:360-dx;
latitude_v  = -90+0.5*dx:dx:90-0.5*dx;

[lat_u,lon_u] = meshgrid(latitude_u,longitude_u);
[lat_v,lon_v] = meshgrid(latitude_v,longitude_v);
[lat_z,lon_z] = meshgrid(latitude_z,longitude_z);

% Convert longitude/latitude from degree to radian
d2r   = pi/180.0;
lon_u = lon_u*d2r;
lat_u = lat_u*d2r;
lon_v = lon_v*d2r;
lat_v = lat_v*d2r;
lon_z = lon_z*d2r;
lat_z = lat_z*d2r;

dlambda = dx*d2r;
dtheta  = dy*d2r;

% compute cos(theta)
cosLatU = cos(lat_u);
cosLatV = cos(lat_v);
cosLatZ = cos(lat_z);

% Reset pole
cosLatZ(:,1  ) = 0.25*cosLatV(:,1);
cosLatZ(:,end) = 0.25*cosLatV(:,end);

% % Plot V grid
% [x,y,z] = sph2cart(lon_v,lat_v,a);
% surf(x,y,z);

% Get grid size
nx_u = size(lon_u,1);
ny_u = size(lon_u,2);
nx_v = size(lon_v,1);
ny_v = size(lon_v,2);
nx_z = size(lon_z,1);
ny_z = size(lon_z,2);

% Initial fields with Rossby-Haurwitz Wave
[u,v,Z] = Haurwitz(a,Omega,g,lon_u,lat_u,lon_v,lat_v,lon_z,lat_z);

% IAP transformation
h                 = sqrt(Z);
hip1(1:nx_z-1,:)  = h(2:nx_z,:);
hip1(nx_z    ,:)  = h(1,:);
hjp1(:,1:ny_z-1)  = h(:,2:ny_z);
hjp1(:,ny_z    )  = 0;

hOnU              = 0.5*(h+hip1); % h on u grid
hOnV_temp         = 0.5*(h+hjp1); % h on v grid
hOnV              = hOnV_temp(:,1:ny_v);

u(:,1  )          = 0;
u(:,end)          = 0;

U                 = hOnU.*u;
V                 = hOnV.*v;

total_energy0     = sum(sum(U.*U.*cosLatU))+sum(sum(V.*V.*cosLatV))+sum(sum(Z.*Z.*cosLatZ));
total_mass0       = sum(sum(Z.*cosLatZ));

% Compute the coefficient for L operator
coefU_x = 0.25./(a*cosLatU*dlambda); %coefU_x = 1/(2*a*cos(theta_u)*2dx)
coefU_y = 0.25./(a*cosLatU*dtheta ); %coefU_y = 1/(2*a*cos(theta_u)*2dy)
coefV_x = 0.25./(a*cosLatV*dlambda); %coefV_x = 1/(2*a*cos(theta_v)*2dx)
coefV_y = 0.25./(a*cosLatV*dtheta ); %coefV_x = 1/(2*a*cos(theta_v)*2dy)
coefZ_x = 0.5 ./(a*cosLatZ*dlambda); %coefZ_x = 1/(  a*cos(theta_z)*2dx)
coefZ_y = 0.5 ./(a*cosLatZ*dtheta ); %coefZ_x = 1/(  a*cos(theta_z)*2dy)

% Output the initial status
int_step_num = ceil(run_time/time_step);
output_count = 0;
output_num   = ceil(run_time/history_interval)+1;

output_netCDF(output_num,output_count,output_precision,...
              U,V,Z,dlambda,dtheta,a,Omega,g,...
              lon_u,lon_v,lon_z,lat_u,lat_v,lat_z,...
              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z)
          
ti_start = clock;
for nt = 1:int_step_num
    % Set previous tend and status for leap-frog
    if nt>1
        U_nm1  = U;
        V_nm1  = V;
        Z_nm1  = Z;
        LU_nm1 = LU;
        LV_nm1 = LV;
        LZ_nm1 = LZ;
        U      = U_np1;
        V      = V_np1;
        Z      = Z_np1;
    else
        U_nm1  = 0;
        V_nm1  = 0;
        Z_nm1  = 0;
        LU_nm1 = 0;
        LV_nm1 = 0;
        LZ_nm1 = 0;
    end
    
    if split_scheme>0
        [tau_n,U_np1,V_np1,Z_np1,LU,LV,LZ] = split_integrator(split_scheme,split_num,IntSch,nt,fast_dt,slow_dt,...
                                                              U_nm1,V_nm1,Z_nm1,LU_nm1,LV_nm1,LZ_nm1,...
                                                              U,V,Z,dlambda,dtheta,a,Omega,g,lat_u,lat_v,lat_z,...
                                                              nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                              coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
    else
        pass = 0;
        [tau_n,U_np1,V_np1,Z_np1,LU,LV,LZ] = integrator(pass,time_step,IntSch,nt,...
                                                        U_nm1,V_nm1,Z_nm1,LU_nm1,LV_nm1,LZ_nm1,...
                                                        U,V,Z,dlambda,dtheta,a,Omega,g,...
                                                        lat_u,lat_v,lat_z,...
                                                        nx_u,ny_u,nx_v,ny_v,nx_z,ny_z,...
                                                        coefU_x,coefU_y,coefV_x,coefV_y,coefZ_x,coefZ_y);
    end
    
    integral_time = nt*time_step;
    
    % Output
	if rem(integral_time,history_interval)==0 && integral_time>=history_interval
        ti_end = clock;
        
        output_count = output_count+1;
        output_netCDF(output_num,output_count,output_precision,...
                      U,V,Z,dlambda,dtheta,a,Omega,g,...
                      lon_u,lon_v,lon_z,lat_u,lat_v,lat_z,...
                      nx_u,ny_u,nx_v,ny_v,nx_z,ny_z)
                  
        total_energy       = sum(sum(U.*U.*cosLatU))+sum(sum(V.*V.*cosLatV))+sum(sum(Z.*Z.*cosLatZ));
        total_mass         = sum(sum(Z.*cosLatZ));
        energy_change_rate = (total_energy-total_energy0)/total_energy0; %ECR
        mass_change_rate   = (total_mass-total_mass0)/total_mass0;       %MCR
        
        disp(['Output        = ',num2str(integral_time/history_interval),'/',num2str(output_num-1)]);
        disp(['tau_n         = ',num2str(tau_n)])
        disp(['Total Energy  = ',num2str(total_energy)]);
        disp(['Total Mass    = ',num2str(total_mass)]);
        disp(['ECR           = ',num2str(energy_change_rate)]);
        disp(['MCR           = ',num2str(mass_change_rate)]);
        disp(['Integral time = ',num2str(etime(ti_end,ti_start))]);
        disp('                                     ');
        
        ti_start = clock;
	end
end

time_end = clock;
disp(['It took ',num2str(etime(time_end,time_start)),' seconds to run SBPEM'])