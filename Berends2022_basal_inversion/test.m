clc
clear all
close all

time   = ncread('exp_I_target_40km/help_fields_ANT.nc','time'); ti = length(time);
phi_ex = ncread('exp_I_target_40km/help_fields_ANT.nc','phi_fric',[1,1,ti],[Inf,Inf,1]);
Hi_ex  = ncread('exp_I_target_40km/restart_ANT.nc'    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);

foldernames = {...
  '../results_20220503_001',...
  '../results_20220503_002',...
  '../results_20220503_003',...
  '../results_20220503_004',...
  '../results_20220503_005',...
  '../results_20220503_006'};

for fi = 1:length(foldernames)
  
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  
  if ~exist(filename,'file'); continue; end

  results(fi).time     = ncread(filename,'time');
  results(fi).RMSE_phi = zeros(size(results(fi).time));

  for ti = 1:length(results(fi).time)
    phi = ncread(filename,'phi_fric',[1,1,ti],[Inf,Inf,1]);
    dphi = log( phi ./ phi_ex);
    dphi(Hi_ex==0) = 0;
    results(fi).RMSE_phi( ti) = sqrt( mean( dphi(:).^2 ));
  end

  plot(results(fi).time,results(fi).RMSE_phi); hold on
  
end

legend(foldernames)

figure; imagesc3(dphi);