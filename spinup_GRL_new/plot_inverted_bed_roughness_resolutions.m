clc
clear all
close all

foldernames = {...
  'phase1_calibration/hybrid_ZoetIverson_40km',...          %  1
  'phase1_calibration/hybrid_ZoetIverson_30km',...          %  2
  'phase1_calibration/hybrid_ZoetIverson_20km'};            %  5

% foldernames = {...
%   'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_40km',...          %  1
%   'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_30km',...          %  2
%   'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_20km',...          %  3
%   'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_16km',...          %  4
%   'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_10km'};            %  5
% 
% foldernames = {...
%   'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_40km',...          %  1
%   'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_30km',...          %  2
%   'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_20km',...          %  3
%   'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_16km',...          %  4
%   'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_10km'};            %  5

cmap_phi  = parula(256);
cmap_Hs   = itmap(16);

cmap_dphi = jet(256);
cmap_dHs  = redbluemap( 32);

cmap_u    = sentinelmap( 256);
cmap_du   = jet( 32);

clim_phi  = [0.1,30];
clim_Hs   = [0,3200];

clim_dphi = [0.1,10];
clim_dHs  = [-250,250];

clim_u    = [1.0,1000];
clim_du   = [0.1,10];

%% Read data
  
for fi = 1:length(foldernames)
  
  filename_restart     = [foldernames{fi} '/restart_GRL.nc'];
  filename_help_fields = [foldernames{fi} '/help_fields_GRL.nc'];
  
  % Read model output
  results(fi).x         = ncread( filename_restart    ,'x');
  results(fi).y         = ncread( filename_restart    ,'y');
  results(fi).time      = ncread( filename_restart    ,'time'); ti = length(results(fi).time);
  results(fi).Hi        = ncread( filename_restart    ,'Hi'       ,[1,1,ti],[Inf,Inf,1]);
  results(fi).Hb        = ncread( filename_restart    ,'Hb'       ,[1,1,ti],[Inf,Inf,1]);
  results(fi).Hs        = ncread( filename_restart    ,'Hs'       ,[1,1,ti],[Inf,Inf,1]);
  results(fi).phi_fric  = ncread( filename_help_fields,'phi_fric' ,[1,1,ti],[Inf,Inf,1]);
  results(fi).uabs_surf = ncread( filename_help_fields,'uabs_surf',[1,1,ti],[Inf,Inf,1]);
  
  % Thickness above flotation
  ice_density           =  910.0;
  seawater_density      = 1028.0;
  results(fi).TAF       = results(fi).Hi - max(0, (-results(fi).Hb) * (seawater_density / ice_density));
end

%% Observed geometry & velocity

% Geometry
for fi = 1: length( foldernames)
  
  filename_restart     = [foldernames{fi} '/restart_GRL.nc'];
  filename_restart = strrep( filename_restart,'phase1a_finecalibration','phase1_calibration');
  filename_restart = strrep( filename_restart,'phase4_holocene','phase1_calibration');
  filename_restart = strrep( filename_restart,'phase5_historicalperiod','phase1_calibration');
  filename_restart = strrep( filename_restart,'PMIP3ens_','');
  filename_restart = strrep( filename_restart,'HadCM3_','');
  filename_restart = strrep( filename_restart,'CCSM_','');
  
  results(fi).Hs_target = ncread( filename_restart,'Hs',[1,1,1],[Inf,Inf,1]);
  results(fi).dHs = results(fi).Hs_target - results(fi).Hs;
end

% Velocity
filename = '/Users/berends/Documents/Datasets/Greenland_velocity/Greenland_ice_velocities_5km.nc';
x = ncread( filename,'x');
y = ncread( filename,'y');
[xg,yg] = meshgrid( y,x);
u = ncread( filename,'u_surf');
v = ncread( filename,'v_surf');
uabs = sqrt( u.^2 + v.^2);
for fi = 1: length( foldernames)
  [xg2,yg2] = meshgrid( results( fi).y, results( fi).x);
  results( fi).uabs_surf_target = interp2( xg, yg, uabs, xg2, yg2);
  results( fi).duabs_surf       = results( fi).uabs_surf ./ results( fi).uabs_surf_target;
end

%% Set up GUI

wa = 150;
ha = size(results(1).Hi,2) * wa / size(results(1).Hi,1);

margins_hor = [110,25,5,5,5,5,110];
margins_ver = [25,5,5,50];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[100,100,wf,hf],'color','w');

H.Ax = zeros(nay,nax);

for i = 1:nax
  for j = 1:nay
    
    x = sum(margins_hor(1:i)) + (i-1)*wa;
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*ha;

    ax = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24,...
      'xtick',[],'ytick',[],'xaxislocation','top',...
      'xlim',[min(results(1).x),max(results(1).x)] * 0.95,...
      'ylim',[min(results(1).y),max(results(1).y)] * 0.95);

    H.Ax(j,i) = ax;

  end
end

for j = 1:nax
  colormap(H.Ax(1,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi,'colorscale','log');

  colormap(H.Ax(2,j),cmap_dHs);
  set(H.Ax(2,j),'clim',clim_dHs);

  colormap(H.Ax(3,j),cmap_du);
  set(H.Ax(3,j),'clim',clim_du,'colorscale','log');
%   colormap(H.Ax(3,j),cmap_u);
%   set(H.Ax(3,j),'clim',clim_u,'colorscale','log');
end

xlabel( H.Ax(1,1),'Observed');
xlabel( H.Ax(1,2),'40 km');
xlabel( H.Ax(1,3),'30 km');
xlabel( H.Ax(1,4),'20 km');
xlabel( H.Ax(1,5),'16 km');
xlabel( H.Ax(1,6),'10 km');

%% Colorbars
pos = get(H.Ax(1,nax),'position');
H.Cbar1 = colorbar( H.Ax(1,nax),'location','eastoutside','ticks',[0.1,0.3,1,3,10,30]);
set(H.Ax(1,nax),'position',pos);
ylabel(H.Cbar1,['Till friction angle (' char(176) ')']);

pos = get(H.Ax(2,nax),'position');
H.Cbar2 = colorbar( H.Ax(2,nax),'location','eastoutside');
set(H.Ax(2,nax),'position',pos);
ylabel(H.Cbar2,'\Delta surface elevation (m)');

pos = get(H.Ax(3,nax),'position');
H.Cbar3 = colorbar( H.Ax(3,nax),'location','eastoutside');
set(H.Ax(3,nax),'position',pos);
ylabel(H.Cbar3,'Surface velocity error (ratio)');

%% Bar between observations and model results

dyp = 15;
axbar = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1)+wa+margins_hor(2)/2, dyp, 10, hf-2*dyp],...
  'xtick',[],'ytick',[],'xlim',[0,1],'ylim',[0,1]);
axbar.XAxis.Visible = 'off';
axbar.YAxis.Visible = 'off';
line('parent',axbar,'xdata',[0.5,0.5],'ydata',[0,1],'linestyle','--','linewidth',2);

%% Plot osberved geometry & velocity

% Colormaps
colormap(H.Ax(2,1),cmap_Hs);
set(H.Ax(2,1),'clim',clim_Hs);

colormap(H.Ax(3,1),cmap_u);
set(H.Ax(3,1),'clim',clim_u,'colorscale','log');

% Colorbars
pos = get(H.Ax(2,1),'position');
H.Cbar4 = colorbar( H.Ax(2,1),'location','westoutside');
set(H.Ax(2,1),'position',pos);
ylabel(H.Cbar4,'Surface elevation (m)');

pos = get(H.Ax(3,1),'position');
H.Cbar5 = colorbar( H.Ax(3,1),'location','westoutside','ticks',[1,10,100,1000]);
set(H.Ax(3,1),'position',pos,'clim',[1,2500]);
ylabel(H.Cbar5,'Surface velocity (m/yr)');

% White patch to hide axeses
xlim = get(H.Ax(1,1),'xlim');
ylim = get(H.Ax(1,1),'ylim');
image('parent',H.Ax( 1,1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 2,1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 3,1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));

% Surface elevation
image('parent',H.Ax( 2,1),'xdata',results( end).x,'ydata',results( end).y,'cdata',results( end).Hs_target','cdatamapping','scaled',...
  'alphadata',double( results( end).TAF' > 0));

% Surface velocity
image('parent',H.Ax( 3,1),'xdata',results( end).x,'ydata',results( end).y,'cdata',results( end).uabs_surf_target','cdatamapping','scaled',...
  'alphadata',double( results( end).TAF' > 0));

%% Plot results

a = 0.5;

% Anomalies
for j = 1: length( foldernames)
  
  % White patch to hide axeses
  xlim = get(H.Ax(1,j+1),'xlim');
  ylim = get(H.Ax(1,j+1),'ylim');
  image('parent',H.Ax( 1,j+1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
  image('parent',H.Ax( 2,j+1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
  image('parent',H.Ax( 3,j+1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
  
  % Till friction angle
  image('parent',H.Ax( 1,j+1),'xdata',results( j).x,'ydata',results( j).y,'cdata',results( j).phi_fric','cdatamapping','scaled',...
    'alphadata',a + (1-a) * double( results( j).TAF' > 0));
  
  % Surface elevation error
  image('parent',H.Ax( 2,j+1),'xdata',results( j).x,'ydata',results( j).y,'cdata',results( j).dHs','cdatamapping','scaled',...
    'alphadata',a + (1-a) * double( results( j).TAF' > 0));
  
  % Surface velocity error
  image('parent',H.Ax( 3,j+1),'xdata',results( j).x,'ydata',results( j).y,'cdata',results( j).duabs_surf','cdatamapping','scaled',...
    'alphadata',a + (1-a) * double( results( j).TAF' > 0));
  
end