clc
clear all
close all

filename_src = 'exp_II_target_2km/help_fields_ANT.nc';
filename_dst = 'exp_II_target_2km/BIV_target_velocity.nc';

%% Create NetCDF template
f = ncinfo( filename_src);

f.Dimensions = f.Dimensions(1:2);

for vi = 1:length(f.Variables)
  if     strcmpi(f.Variables(vi).Name,'u_surf')
    var_u = f.Variables(vi);
  elseif strcmpi(f.Variables(vi).Name,'v_surf')
    var_v = f.Variables(vi);
  end
end

var_u.Dimensions = f.Dimensions;
var_v.Dimensions = f.Dimensions;
var_u.Size = [f.Dimensions(1).Length, f.Dimensions(2).Length];
var_v.Size = [f.Dimensions(1).Length, f.Dimensions(2).Length];

f.Variables = f.Variables(1:2);
f.Variables(end+1) = var_u;
f.Variables(end+1) = var_v;

%% Create file
if exist(filename_dst,'file')
  delete(filename_dst)
end

ncwriteschema( filename_dst,f);

%% Read and write data

x    = ncread( filename_src,'x');
y    = ncread( filename_src,'y');
time = ncread( filename_src,'time'); ti = length(time);
u    = ncread( filename_src,'u_surf',[1,1,ti],[Inf,Inf,1]);
v    = ncread( filename_src,'v_surf',[1,1,ti],[Inf,Inf,1]);

ncwrite( filename_dst,'x'     ,x);
ncwrite( filename_dst,'y'     ,y);
ncwrite( filename_dst,'u_surf',u);
ncwrite( filename_dst,'v_surf',v);