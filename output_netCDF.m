function output_netCDF(MESH,STATE,history_interval,output_count,output_precision,output_num)

f_out   = 'output.nc';

% Inverse IAP transformation
h                     = sqrt(STATE.Z);
hip1(1:MESH.nx_z-1,:) = h(2:MESH.nx_z,:);
hip1(MESH.nx_z    ,:) = h(1,:);
hjp1(:,1:MESH.ny_z-1) = h(:,2:MESH.ny_z);
hjp1(:,MESH.ny_z    ) = 0;

hOnU      = 0.5*(h+hip1); % h on u grid
hOnV_temp = 0.5*(h+hjp1); % h on v grid
hOnV      = hOnV_temp(:,1:MESH.ny_v);

u         = STATE.U./hOnU;
v         = STATE.V./hOnV;

% Write Data into netCDF
west_east        = MESH.nx_u;
south_north      = MESH.ny_u;
south_north_stag = MESH.ny_v;

if output_count==0
    mode           = netcdf.getConstant('NETCDF4');
    mode           = bitor(mode,netcdf.getConstant('CLOBBER'));
    ncid           = netcdf.create(f_out,mode);
    disp(['ncid = ',num2str(ncid)])
    
    % Define Dimensions
    time_dimID             = netcdf.defDim(ncid,'time'            ,output_num);
    west_east_dimID        = netcdf.defDim(ncid,'west_east'       ,west_east);
    south_north_dimID      = netcdf.defDim(ncid,'south_north'     ,south_north);
    south_north_stag_dimID = netcdf.defDim(ncid,'south_north_stag',south_north_stag);
    
    % Define Attribute
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model_Source'     ,'SBPEM written by Zhou Lilong')
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dlambda'          ,MESH.dlambda)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dtheta'           ,MESH.dtheta)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'earth_radius'     ,MESH.a)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Omega'            ,MESH.Omega)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'g'                ,MESH.g)
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history_interval' ,history_interval)
    
    % Define Variables
    XLONG_U_id = netcdf.defVar(ncid,'XLONG_U',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLONG_U_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_U_id,'description','longitude for u');
    
    XLAT_U_id = netcdf.defVar(ncid,'XLAT_U',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLAT_U_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_U_id,'description','latitude for u');
    
    XLONG_V_id = netcdf.defVar(ncid,'XLONG_V',output_precision,[west_east_dimID,south_north_stag_dimID]);
    netcdf.putAtt(ncid,XLONG_V_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_V_id,'description','longitude for v');
    
    XLAT_V_id = netcdf.defVar(ncid,'XLAT_V',output_precision,[west_east_dimID,south_north_stag_dimID]);
    netcdf.putAtt(ncid,XLAT_V_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_V_id,'description','latitude for v');
    
    XLONG_M_id = netcdf.defVar(ncid,'XLONG_M',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLONG_M_id,'units','degrees longitude');
    netcdf.putAtt(ncid,XLONG_M_id,'description','longitude for z');
    
    XLAT_M_id = netcdf.defVar(ncid,'XLAT_M',output_precision,[west_east_dimID,south_north_dimID]);
    netcdf.putAtt(ncid,XLAT_M_id,'units','degrees latitude');
    netcdf.putAtt(ncid,XLAT_M_id,'description','latitude for z');
    
    U_id = netcdf.defVar(ncid,'U',output_precision,[west_east_dimID,south_north_dimID,time_dimID]);
    netcdf.putAtt(ncid,U_id,'units','m/s');
    netcdf.putAtt(ncid,U_id,'description','u wind component');
    
    V_id = netcdf.defVar(ncid,'V',output_precision,[west_east_dimID,south_north_stag_dimID,time_dimID]);
    netcdf.putAtt(ncid,V_id,'units','m/s');
    netcdf.putAtt(ncid,V_id,'description','v wind component');
    
    Z_id = netcdf.defVar(ncid,'Z',output_precision,[west_east_dimID,south_north_dimID,time_dimID]);
    netcdf.putAtt(ncid,Z_id,'units','m^2/s^2');
    netcdf.putAtt(ncid,Z_id,'description','geopotential height');
    
    % Put Variables
    netcdf.putVar(ncid, XLONG_U_id ,MESH.lon_u);
    netcdf.putVar(ncid, XLAT_U_id  ,MESH.lat_u);
    netcdf.putVar(ncid, XLONG_V_id ,MESH.lon_v);
    netcdf.putVar(ncid, XLAT_V_id  ,MESH.lat_v);
    netcdf.putVar(ncid, XLONG_M_id ,MESH.lon_z);
    netcdf.putVar(ncid, XLAT_M_id  ,MESH.lat_z);
    netcdf.putVar(ncid, U_id       ,[0,0,0],[west_east,south_north     ,1],u      );
    netcdf.putVar(ncid, V_id       ,[0,0,0],[west_east,south_north_stag,1],v      );
    netcdf.putVar(ncid, Z_id       ,[0,0,0],[west_east,south_north     ,1],STATE.Z);
    
    netcdf.close(ncid)
    
else
    ncid = netcdf.open(f_out,'WRITE');
    
    U_id = netcdf.inqVarID(ncid,'U');
    V_id = netcdf.inqVarID(ncid,'V');
    Z_id = netcdf.inqVarID(ncid,'Z');
    
    netcdf.reDef(ncid)
    
    netcdf.putVar(ncid, U_id     ,[0,0,output_count],[west_east,south_north     ,1],u      );
    netcdf.putVar(ncid, V_id     ,[0,0,output_count],[west_east,south_north_stag,1],v      );
    netcdf.putVar(ncid, Z_id     ,[0,0,output_count],[west_east,south_north     ,1],STATE.Z);
    
    netcdf.close(ncid)
end

