clc
clear all
close all

filename = 'BIVMIP_B_inv_10km_perfect/restart_ANT.nc';

time = ncread( filename,'time');

Hi_ref = ncread( filename,'Hi',[1,1,1],[Inf,Inf,1]);

RMSE = [];

for ti = 1: length( time)
  Hi = ncread( filename,'Hi',[1,1,ti],[Inf,Inf,1]);
  RMSE( end+1) = sqrt( mean( (Hi(:) - Hi_ref(:)).^2 ));
end

plot( RMSE)