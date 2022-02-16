clc
clear all
close all

filename1 = 'BIVMIP_C_inv_5km_perfect/inverted_bed_roughness.nc';
filename2 = 'BIVMIP_C_inv_5km_SMB_lo/inverted_bed_roughness.nc';

if exist( filename2,'file')
  delete( filename2)
end

copyfile( filename1, filename2)

phi = ncread('BIVMIP_C_inv_5km_SMB_lo/help_fields_ANT.nc','phi_fric');
phi = phi(:,:,end);

ncwrite(filename2,'phi_fric',phi);