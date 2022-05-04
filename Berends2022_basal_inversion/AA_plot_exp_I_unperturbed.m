clc
clear all
close all

foldernames = {...
  'exp_I_target_40km',...
  'exp_I_target_20km',...
  'exp_I_target_10km',...
  'exp_I_inv_40km_unperturbed',...
  'exp_I_inv_20km_unperturbed',...
  'exp_I_inv_10km_unperturbed'};

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
end

for fi = 4:6
  results(fi).dHs       = results(fi).Hs       - results(fi-3).Hs;
  results(fi).dphi_fric = results(fi).phi_fric ./ results(fi-3).phi_fric;
  results(fi).du        = results(fi).uabs     - results(fi-3).uabs;
end

%% Set up GUI

wa = 250;
ha = 250;

margins_hor = [120,50,5,5,120];
margins_ver = [25,25,25,50];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[100,100,wf,hf],'color','w');

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

for j = 1:1
  colormap(H.Ax(1,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi);
  
  colormap(H.Ax(2,j),cmap_Hs);
  set(H.Ax(2,j),'clim',clim_Hs);
  
  colormap(H.Ax(3,j),cmap_u);
  set(H.Ax(3,j),'clim',clim_u);
end

for j = 2:4
  colormap(H.Ax(1,j),cmap_dphi);
  set(H.Ax(1,j),'clim',clim_dphi,'colorscale','log');
  
  colormap(H.Ax(2,j),cmap_dHs);
  set(H.Ax(2,j),'clim',clim_dHs);
  
  colormap(H.Ax(3,j),cmap_du);
  set(H.Ax(3,j),'clim',clim_du);
end

xlabel( H.Ax(1,1),'Target (10 km)');
xlabel( H.Ax(1,2),'Unperturbed (40 km)');
xlabel( H.Ax(1,3),'20 km');
xlabel( H.Ax(1,4),'10 km');

for j = 1:4
  
  r = 0.86;
  set(H.Ax(1,j),'xlim',[min(results(1).x),max(results(1).x)]*r,'ylim',[min(results(1).y),max(results(1).y)]*r);
  set(H.Ax(2,j),'xlim',[min(results(1).x),max(results(1).x)]*r,'ylim',[min(results(1).y),max(results(1).y)]*r);
  set(H.Ax(3,j),'xlim',[min(results(1).x),max(results(1).x)]*r,'ylim',[min(results(1).y),max(results(1).y)]*r);
  
end

%% Line between first and second column

x = margins_hor(1) + wa + margins_hor(2)/2 - 10;
pos = get(H.Fig,'position');
ax = axes('parent',H.Fig,'units','pixels','position',[x,25,20,pos(4)-50],'xlim',[0,1],'ylim',[0,1]);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
line('parent',ax,'xdata',[0.5,0.5],'ydata',[0,1],'linestyle','--')

%% Colorbars
pos = get(H.Ax(1,1),'position');
H.Cbarl1 = colorbar(H.Ax(1,1),'location','westoutside');
set(H.Ax(1,1),'position',pos);
ylabel( H.Cbarl1,['Till friction angle (' char(176) ')'])

pos = get(H.Ax(2,1),'position');
H.Cbarl2 = colorbar(H.Ax(2,1),'location','westoutside');
set(H.Ax(2,1),'position',pos);
ylabel( H.Cbarl2,'Surface elevation (m)')

pos = get(H.Ax(3,1),'position');
H.Cbarl3 = colorbar(H.Ax(3,1),'location','westoutside');
set(H.Ax(3,1),'position',pos);
ylabel( H.Cbarl3,'Surface velocity (m yr^{-1})')
set(H.Cbarl3,'ticks',[1,10,100,1000]);

pos = get(H.Ax(1,4),'position');
H.Cbarr1 = colorbar(H.Ax(1,4),'location','eastoutside');
set(H.Ax(1,4),'position',pos);
ylabel( H.Cbarr1,['\Delta till friction angle (' char(176) ')'])

pos = get(H.Ax(2,4),'position');
H.Cbarr2 = colorbar(H.Ax(2,4),'location','eastoutside');
set(H.Ax(2,4),'position',pos);
ylabel( H.Cbarr2,'\Delta surface elevation (m)')

pos = get(H.Ax(3,4),'position');
H.Cbarr3 = colorbar(H.Ax(3,4),'location','eastoutside');
set(H.Ax(3,4),'position',pos);
ylabel( H.Cbarr3,'\Delta surface velocity (m yr^{-1})')
  
%% Plot results - target
x = results(3).x;
y = results(3).y;

% Blank white image to cover the axes lines (because they're ugly)
cdata = zeros(length(y),length(x),3)+1;
image('parent',H.Ax(1,1),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
image('parent',H.Ax(2,1),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
image('parent',H.Ax(3,1),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);

% Top row: till friction angle
cdata = results(3).phi_fric';
adata = zeros(size(cdata));
adata( results(3).Hi'>0) = 1;
image('parent',H.Ax(1,1),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Middle row: surface elevation
cdata = results(3).Hs';
adata = zeros(size(cdata));
adata( results(3).Hi'>0) = 1;
image('parent',H.Ax(2,1),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Bottom row: surface velocity
cdata = results(3).uabs';
adata = zeros(size(cdata));
adata( results(3).Hi'>0) = 1;
image('parent',H.Ax(3,1),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
set(H.Ax(3,1),'colorscale','log');
  
%% Plot results - inverted
for j = 1:3
  
  jr = j+3;
  ja = j+1;
  
  x = results(jr).x;
  y = results(jr).y;

  % Blank white image to cover the axes lines (because they're ugly)
  cdata = zeros(length(y),length(x),3)+1;
  image('parent',H.Ax(1,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  image('parent',H.Ax(2,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  image('parent',H.Ax(3,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);

  % Top row: till friction angle
  cdata = results(jr).dphi_fric';
  adata = zeros(size(cdata));
  adata( results(jr).Hi'>10) = 1;
  image('parent',H.Ax(1,ja),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

  % Middle row: surface elevation
  cdata = results(jr).dHs';
  adata = zeros(size(cdata));
  adata( results(jr).Hi'>10) = 1;
  image('parent',H.Ax(2,ja),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

  % Bottom row: surface velocity
  cdata = results(jr).du';
  adata = zeros(size(cdata));
  adata( results(jr).Hi'>10) = 1;
  image('parent',H.Ax(3,ja),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
end