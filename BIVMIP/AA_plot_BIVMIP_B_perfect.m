clc
clear all
close all

foldernames = {...
  'BIVMIP_B_perfect_10km',...
  'BIVMIP_B_inv_40km_perfect',...
  'BIVMIP_B_inv_40km_perfect',...
  'BIVMIP_B_inv_40km_perfect'};

cmap_phi = parula(256);
cmap_Hs  = itmap(16);

clim_phi = [0,6];
clim_Hs  = [0,2700];

%% Read data
for fi = 1:length(foldernames)
  
  filename_restart     = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields = [foldernames{fi} '/help_fields_ANT.nc'];
  
  x         = ncread( filename_restart    ,'x');
  y         = ncread( filename_restart    ,'y');
  time      = ncread( filename_restart    ,'time'); ti = length(time);
  Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
  Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
  phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  
  results(fi).x        = x;
  results(fi).y        = y;
  results(fi).Hi       = Hi;
  results(fi).Hb       = Hb;
  results(fi).Hs       = Hs;
  results(fi).phi_fric = phi_fric;
end

%% Set up GUI

wa = 250;
ha = 250;

margins_hor = [5,5,5,5,120];
margins_ver = [25,25,50];

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

for j = 1:4
  colormap(H.Ax(1,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi);
  
  colormap(H.Ax(2,j),cmap_Hs);
  set(H.Ax(2,j),'clim',clim_Hs);
end

xlabel( H.Ax(1,1),'Target (10 km)');
xlabel( H.Ax(1,2),'Inverted (40 km)');
xlabel( H.Ax(1,3),'Inverted (20 km)');
xlabel( H.Ax(1,4),'Inverted (10 km)');

pos = get(H.Ax(1,4),'position');
H.Cbar1 = colorbar(H.Ax(1,4),'location','eastoutside');
set(H.Ax(1,4),'position',pos);
ylabel( H.Cbar1,['Till friction angle (' char(176) ')'])

pos = get(H.Ax(2,4),'position');
H.Cbar2 = colorbar(H.Ax(2,4),'location','eastoutside');
set(H.Ax(2,4),'position',pos);
ylabel( H.Cbar2,'Surface elevation (m)')

%% Plot data

for j = 1:4
  
  r = 0.86;
  set(H.Ax(1,j),'xlim',[min(results(1).x),max(results(1).x)]*r,'ylim',[min(results(1).y),max(results(1).y)]*r);
  set(H.Ax(2,j),'xlim',[min(results(1).x),max(results(1).x)]*r,'ylim',[min(results(1).y),max(results(1).y)]*r);
  
  x = results(j).x; dx = x(2) - x(1);
  y = results(j).y;
  
  % Blank white image to cover the axes lines (because they're ugly)
  cdata = zeros(length(y),length(x),3)+1;
  image('parent',H.Ax(1,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  image('parent',H.Ax(2,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  
  % Top row: till friction angle
  cdata = results(j).phi_fric';
  adata = zeros(size(cdata));
  adata( results(j).Hi'>0) = 1;
  image('parent',H.Ax(1,j),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
  
  % Bottom row: surface elevation
  cdata = results(j).Hs';
  adata = zeros(size(cdata));
  adata( results(j).Hi'>0) = 1;
  image('parent',H.Ax(2,j),'xdata',x,'ydata',y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
end