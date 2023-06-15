clc
clear all
close all

% Plot errors in elevation and velocity at the end of the historical period

run_name = 'hybrid_ZoetIverson_PMIP3ens_16km';

cmap_Hs   = itmap( 17);
cmap_dHs  = redbluemap( 17);

cmap_u    = sentinelmap( 17);
cmap_du   = jet( 17);

clim_Hs   = [0,3200];
clim_dHs  = [-250,250];

clim_u    = [1.0,1500];
clim_du   = [0.1,10];

%% Read observed geometry and velocity

% Geometry
filename = ['phase1_calibration/' strrep(run_name,'PMIP3ens_','') '/restart_GRL.nc'];
PD.x  = ncread( filename,'x');
PD.y  = ncread( filename,'y');
PD.Hs = ncread( filename,'Hs',[1,1,1],[Inf,Inf,1]);
PD.Hi = ncread( filename,'Hi',[1,1,1],[Inf,Inf,1]);

% Velocity
filename = '/Users/berends/Documents/Datasets/Greenland_velocity/Greenland_ice_velocities_5km.nc';
x = ncread( filename,'x');
y = ncread( filename,'y');
[xg,yg] = meshgrid( y,x);
u = ncread( filename,'u_surf');
v = ncread( filename,'v_surf');
uabs = sqrt( u.^2 + v.^2);
[PD.xg,PD.yg] = meshgrid( PD.y, PD.x);
PD.u = interp2( xg,yg,uabs,PD.xg,PD.yg);

%% Read modelled geometry and velocity

% Geometry
filename = ['phase4_holocene/' run_name '/restart_GRL.nc'];
% filename = ['phase5_historicalperiod/' run_name '/restart_GRL.nc'];
time    = ncread( filename,'time'); ti = length( time);
spin.Hi = ncread( filename,'Hi',[1,1,ti],[Inf,Inf,1]);
spin.Hs = ncread( filename,'Hs',[1,1,ti],[Inf,Inf,1]);

% Velocity
filename = ['phase4_holocene/' run_name '/help_fields_GRL.nc'];
% filename = ['phase5_historicalperiod/' run_name '/help_fields_GRL.nc'];
time    = ncread( filename,'time'); ti = length( time);
spin.u  = ncread( filename,'uabs_surf',[1,1,ti],[Inf,Inf,1]);

%% Remove Ellesmere Island
grid.x = ncread( filename,'x'); grid.nx = length( grid.x);
grid.y = ncread( filename,'y'); grid.ny = length( grid.y);
grid.lambda_M    = -45.0;
grid.phi_M       = 90.0;
grid.beta_stereo = 70.0;

mask_noice = initialise_mask_noice_GRL_remove_Ellesmere( grid);
PD.Hi(   mask_noice == 1) = 0;
PD.Hs(   mask_noice == 1) = 0;
PD.u(    mask_noice == 1) = 0;
spin.Hi( mask_noice == 1) = 0;
spin.Hs( mask_noice == 1) = 0;
spin.u(  mask_noice == 1) = 0;

%% Set up GUI

wa = 200;
ha = size( PD.Hs,2) * wa / size( PD.Hs,1);

margins_hor = [130,15,15,130];
margins_ver = [25,25,50];

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

    ax = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',28,...
      'xtick',[],'ytick',[],'xlim',[0,1],'ylim',[0,1],'xaxislocation','top');

    H.Ax(j,i) = ax;

  end
end

xlabel(H.Ax(1,1),'observed')
xlabel(H.Ax(1,2),'spun-up')
xlabel(H.Ax(1,3),'difference')

% Colormaps
colormap(H.Ax(1,1),cmap_Hs);
colormap(H.Ax(1,2),cmap_Hs);
colormap(H.Ax(1,3),cmap_dHs);
set(H.Ax(1,1),'clim',clim_Hs);
set(H.Ax(1,2),'clim',clim_Hs);
set(H.Ax(1,3),'clim',clim_dHs);
colormap(H.Ax(2,1),cmap_u);
colormap(H.Ax(2,2),cmap_u);
colormap(H.Ax(2,3),cmap_du);
set(H.Ax(2,1),'clim',clim_u,'colorscale','log');
set(H.Ax(2,2),'clim',clim_u,'colorscale','log');
set(H.Ax(2,3),'clim',clim_du,'colorscale','log');

% Colorbars
pos = get(H.Ax(1,1),'position');
H.Cbar1 = colorbar(H.Ax(1,1),'location','westoutside');
set(H.Ax(1,1),'position',pos);
ylabel(H.Cbar1,'Surface elevation (m)')

pos = get(H.Ax(1,3),'position');
H.Cbar2 = colorbar(H.Ax(1,3),'location','eastoutside');
set(H.Ax(1,3),'position',pos);
ylabel(H.Cbar2,'\Delta surface elevation (m)')

pos = get(H.Ax(2,1),'position');
H.Cbar3 = colorbar(H.Ax(2,1),'location','westoutside','ticks',[0.1,1,10,100,1000]);
set(H.Ax(2,1),'position',pos);
ylabel(H.Cbar3,'Surface velocity (m/yr)')

pos = get(H.Ax(2,3),'position');
H.Cbar4 = colorbar(H.Ax(2,3),'location','eastoutside');
set(H.Ax(2,3),'position',pos);
ylabel(H.Cbar4,'\Delta surface velocity')

%% Plot results

% White patch to hide axeses
xlim = get(H.Ax(1,1),'xlim');
ylim = get(H.Ax(1,1),'ylim');
image('parent',H.Ax( 1,1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 1,2),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 1,3),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 2,1),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 2,2),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));
image('parent',H.Ax( 2,3),'xdata',xlim*1.1,'ydata',ylim*1.1,'cdata',ones(2,2,3));

% Data
image('parent',H.Ax(1,1),'xdata',[0,1],'ydata',[0,1],'cdata',PD.Hs'           ,'cdatamapping','scaled');
image('parent',H.Ax(1,2),'xdata',[0,1],'ydata',[0,1],'cdata',spin.Hs'         ,'cdatamapping','scaled');
image('parent',H.Ax(1,3),'xdata',[0,1],'ydata',[0,1],'cdata',spin.Hs' - PD.Hs','cdatamapping','scaled');

image('parent',H.Ax(2,1),'xdata',[0,1],'ydata',[0,1],'cdata',PD.u'            ,'cdatamapping','scaled','alphadata',double(spin.Hi'>50));
image('parent',H.Ax(2,2),'xdata',[0,1],'ydata',[0,1],'cdata',spin.u'          ,'cdatamapping','scaled','alphadata',double(spin.Hi'>50));
image('parent',H.Ax(2,3),'xdata',[0,1],'ydata',[0,1],'cdata',spin.u'  ./ PD.u','cdatamapping','scaled','alphadata',double(spin.Hi'>50));

function mask_noice = initialise_mask_noice_GRL_remove_Ellesmere( grid)
  % Prevent ice growth in the Ellesmere Island part of the Greenland domain
  
  mask_noice = zeros( grid.ny, grid.nx);

  % The two endpoints in lat,lon
  pa_latlon = [76.74, -74.79];
  pb_latlon = [82.19, -60.00];

  % The two endpoints in x,y
  [xa,ya] = oblique_sg_projection( pa_latlon(2), pa_latlon(1), grid.lambda_M, grid.phi_M, grid.beta_stereo);
  [xb,yb] = oblique_sg_projection( pb_latlon(2), pb_latlon(1), grid.lambda_M, grid.phi_M, grid.beta_stereo);

  pa = [xa,ya];
  pb = [xb,yb];

  for i = 1: grid.nx
    yl_ab = pa(2) + (grid.x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1));
    for j = 1: grid.ny
      if (grid.y(j) > pa(2) && grid.y(j) > yl_ab && grid.x(i) < pb(1))
        mask_noice( j,i) = 1;
      end
    end
  end
  
  mask_noice = mask_noice';
  
  function [x_IM_P_prime, y_IM_P_prime] = oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, beta_deg)
    
    earth_radius                     = 6.371221E6;
    
    % Convert beta to alpha
    alpha_deg = 90.0 - beta_deg;

    % Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180.0) * phi;
    lambda_P = (pi / 180.0) * lambda;

    % Convert projection parameters to radians:
    lambda_M = (pi / 180.0) * lambda_M_deg;
    phi_M    = (pi / 180.0) * phi_M_deg;
    alpha    = (pi / 180.0) * alpha_deg;

    % See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1.0 + cos(alpha)) / (1.0 + cos(phi_P) * cos(phi_M) * cos(lambda_P - lambda_M) + sin(phi_P) * sin(phi_M));

    % See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  earth_radius * (cos(phi_P) * sin(lambda_P - lambda_M)) * t_P_prime;
    y_IM_P_prime =  earth_radius * (sin(phi_P) * cos(phi_M) - (cos(phi_P) * sin(phi_M)) * cos(lambda_P - lambda_M)) * t_P_prime;

  end
  
end