clc
clear all
close all

help_fields_filename = 'BIVMIP_B_perfect_10km/help_fields_ANT.nc';
BIV_target_filename  = 'BIVMIP_B_velocity_10km.nc';

%% Set up NetCDF template
f = ncinfo( help_fields_filename);

% Dimensions: only x,y
f.Dimensions = f.Dimensions(1:2);

% Variables: x,y,u_surf,v_surf
f.Variables = f.Variables([1,2,7,7]);

% u_surf
f.Variables(3).Name = 'u_surf';
f.Variables(3).Dimensions = f.Dimensions;
f.Variables(3).Size = [f.Dimensions(1).Length, f.Dimensions(2).Length];
f.Variables(3).Attributes(1).Name  = 'long_name';
f.Variables(3).Attributes(1).Value = 'surface velocity in x-direction';
f.Variables(3).Attributes(2).Name  = 'units';
f.Variables(3).Attributes(2).Value = 'm/yr';

% v_surf
f.Variables(4) = f.Variables(3);
f.Variables(4).Name = 'v_surf';

%% Create file, read+write data
if exist(BIV_target_filename,'file')
  delete(BIV_target_filename)
end

ncwriteschema( BIV_target_filename, f);

x      = ncread( help_fields_filename,'x');
y      = ncread( help_fields_filename,'y');
time   = ncread( help_fields_filename,'time');
ti = length(time);

f = ncinfo( help_fields_filename);
u = [];
v = [];
for vi = 1: length(f.Variables)
  if     strcmpi( f.Variables(vi).Name,'u_surf')
    u = ncread( help_fields_filename,'u_surf',[1,1,ti],[Inf,Inf,1]);
  elseif strcmpi( f.Variables(vi).Name,'u_3D')
    u = ncread( help_fields_filename,'u_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  elseif strcmpi( f.Variables(vi).Name,'v_surf')
    v = ncread( help_fields_filename,'v_surf',[1,1,ti],[Inf,Inf,1]);
  elseif strcmpi( f.Variables(vi).Name,'v_3D')
    v = ncread( help_fields_filename,'v_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  end
end

ncwrite( BIV_target_filename,'x'     ,x);
ncwrite( BIV_target_filename,'y'     ,y);
ncwrite( BIV_target_filename,'u_surf',u);
ncwrite( BIV_target_filename,'v_surf',v);