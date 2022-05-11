clc
clear all
close all

foldername_target = 'exp_II_target_5km';

foldernames = {...
  'exp_II_inv_5km_noneq'};            %  8

cmap_phi  = parula(256);
cmap_Hs   = itmap(16);

cmap_dphi = jet(32);
cmap_dHs  = flipud(lbmap(32,'redblue'));

cmap_u    = sentinelmap( 32);
cmap_du   = jet( 32);

clim_phi  = [0,6];
clim_Hs   = [0,2700];

clim_dphi = [0.1,10];
clim_dHs  = [-250,250];

clim_u    = [1.0,1000];
clim_du   = [-250,250];

%% Read data
  
filename_restart     = [foldername_target '/restart_ANT.nc'];
filename_help_fields = [foldername_target '/help_fields_ANT.nc'];

target.x         = ncread( filename_restart    ,'x');
target.y         = ncread( filename_restart    ,'y');
target.time      = ncread( filename_restart    ,'time'); ti = find(target.time==10000);
target.Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
target.Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
target.Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
target.phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
target.u         = ncread( filename_help_fields,'u_surf'  ,[1,1,ti],[Inf,Inf,1]);
target.v         = ncread( filename_help_fields,'v_surf'  ,[1,1,ti],[Inf,Inf,1]);
target.uabs      = sqrt( target.u.^2 + target.v.^2);
  
% Remove the weird artefact in the northwest corner
target.Hi(1,:) = min(target.Hi(1,:));
target.Hs(1,:) = min(target.Hs(1,:));

ice_density      =  910.0;
seawater_density = 1028.0;
target.TAF       = target.Hi - max(0, (-target.Hb) * (seawater_density / ice_density));
  
for fi = 1:length(foldernames)
  
  filename_restart     = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields = [foldernames{fi} '/help_fields_ANT.nc'];
  
  results(fi).x         = ncread( filename_restart    ,'x');
  results(fi).y         = ncread( filename_restart    ,'y');
  results(fi).time      = ncread( filename_restart    ,'time'); ti = length(results(fi).time);
  results(fi).Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  results(fi).u         = ncread( filename_help_fields,'u_surf'  ,[1,1,ti],[Inf,Inf,1]);
  results(fi).v         = ncread( filename_help_fields,'v_surf'  ,[1,1,ti],[Inf,Inf,1]);
  results(fi).uabs      = sqrt( results(fi).u.^2 + results(fi).v.^2);
  
  % Remove the weird artefact in the northwest corner
  results(fi).Hi(1,:) = min(results(fi).Hi(1,:));
  results(fi).Hs(1,:) = min(results(fi).Hs(1,:));
  
  results(fi).dphi_fric = results(fi).phi_fric ./ target.phi_fric;
  results(fi).dHs       = results(fi).Hs       -  target.Hs;
  results(fi).du        = results(fi).uabs     -  target.uabs;
  
  ice_density           =  910.0;
  seawater_density      = 1028.0;
  results(fi).TAF       = results(fi).Hi - max(0, (-results(fi).Hb) * (seawater_density / ice_density));
end

%% Set up GUI

wa = 400;
ha = 100;

margins_hor = [25,25];
margins_ver = [90,90,90,15];

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
      'xtick',[],'ytick',[],'xaxislocation','top');

    H.Ax(j,i) = ax;

  end
end

colormap(H.Ax(1,1),cmap_dphi);
set(H.Ax(1,1),'clim',clim_dphi,'colorscale','log');
pos = get(H.Ax(1,1),'position');
H.Cbar1 = colorbar(H.Ax(1,1),'location','southoutside');
set(H.Ax(1,1),'position',pos);
ylabel(H.Cbar1,['\Delta till friction angle (' char(176) ')']);

colormap(H.Ax(2,1),cmap_dHs);
set(H.Ax(2,1),'clim',clim_dHs);
pos = get(H.Ax(2,1),'position');
H.Cbar2 = colorbar(H.Ax(2,1),'location','southoutside');
set(H.Ax(2,1),'position',pos);
ylabel(H.Cbar2,'\Delta surface elevation (m)');

colormap(H.Ax(3,1),cmap_du);
set(H.Ax(3,1),'clim',clim_du);
pos = get(H.Ax(3,1),'position');
H.Cbar3 = colorbar(H.Ax(3,1),'location','southoutside');
set(H.Ax(3,1),'position',pos);
ylabel(H.Cbar3,'\Delta surface velocity (m yr^{-1})');

%% Plot results

for i = 1:3
  for j = 1:1
    
    x = results(1).x;
    y = results(1).y;
    set(H.Ax(i,j),'xlim',[min(x),max(x)*0.5],'ylim',[min(y),max(y)]);

    % Blank white image to cover the axes lines (because they're ugly)
    cdata = zeros(length(y),length(x),3)+1;
    image('parent',H.Ax(i,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);

  end
end

% till friction angle

ax = H.Ax(1,1);
R  = results(1);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.TAF'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% surface elevation

ax = H.Ax(2,1);
R  = results(1);
cdata = R.dHs';
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

% surface velocity

ax = H.Ax(3,1);
R  = results(1);
cdata = R.du';
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

%% Grounding lines
  
% Construct target grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL_target = contour('parent',H.tempax,'xdata',target.y,'ydata',target.x,'zdata',target.TAF,'levellist',0);
close(H.tempfig);

% Viscosity
ax = H.Ax(2,1);
R  = results(1);
  
% Construct actual grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL = contour('parent',H.tempax,'xdata',R.y,'ydata',R.x,'zdata',R.TAF,'levellist',0);
close(H.tempfig);

% Plot target GL contour
C = C_GL_target;
while ~isempty(C)
  n  = C(2,1);
  Ct = C(:,2:2+n-1);
  C = C(:,2+n:end);
  line('parent',ax,'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
end

% Plot inverted GL contour
C = C_GL;
while ~isempty(C)
  n  = C(2,1);
  Ct = C(:,2:2+n-1);
  C = C(:,2+n:end);
  line('parent',ax,'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','k','linestyle','--');
end
