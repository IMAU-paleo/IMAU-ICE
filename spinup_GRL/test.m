clc
clear all
close all

filename1 = 'phase1_calibration/hybrid_ZoetIverson_20km/help_fields_GRL.nc';
filename2 = 'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_20km/help_fields_GRL.nc';

Hs1 = ncread( filename1,'Hs'      ,[1,1,1],[Inf,Inf,1]);
T1  = ncread( filename1,'T2m_year',[1,1,1],[Inf,Inf,1]);
T1  = T1 + Hs1 * 0.008;

time = ncread( filename2,'time');
nt = find( time==1980);
Hs2 = mean( ncread( filename2,'Hs'      ,[1,1,1],[Inf,Inf,nt]), 3);
T2  = mean( ncread( filename2,'T2m_year',[1,1,1],[Inf,Inf,nt]), 3);
T2  = T2 + Hs2 * 0.008;

figure; imagesc3( T1); set(gca,'clim',[255,280]);
figure; imagesc3( T2); set(gca,'clim',[255,280]);
figure; imagesc3( T2 - T1); set(gca,'clim',[-0.05,0.05]);