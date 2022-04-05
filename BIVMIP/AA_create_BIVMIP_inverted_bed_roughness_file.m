clc
clear all
close all

filename_src = 'BIVMIP_B_perfect_40km/help_fields_ANT.nc';
filename_dst = 'BIVMIP_B_perfect_40km/inverted_bed_roughness.nc';

%% Create NetCDF template
f = ncinfo( filename_src);

f.Dimensions = f.Dimensions(1:2);

for vi = 1:length(f.Variables)
  if     strcmpi(f.Variables(vi).Name,'phi_fric')
    var_phi = f.Variables(vi);
  end
end

var_u.Dimensions = f.Dimensions;
var_u.Size = [f.Dimensions(1).Length, f.Dimensions(2).Length];

f.Variables = f.Variables(1:2);
f.Variables(end+1) = var_u;

%% Create file
if exist(filename_dst,'file')
  delete(filename_dst)
end

ncwriteschema( filename_dst,f);

%% Read and write data

x    = ncread( filename_src,'x');
y    = ncread( filename_src,'y');
time = ncread( filename_src,'time'); ti = length(time);
phi  = ncread( filename_src,'phi_fric',[1,1,ti],[Inf,Inf,1]);

ncwrite( filename_dst,'x'       ,x  );
ncwrite( filename_dst,'y'       ,y  );
ncwrite( filename_dst,'phi_fric',phi);