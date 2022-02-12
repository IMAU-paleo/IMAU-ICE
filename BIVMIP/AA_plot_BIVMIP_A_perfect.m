clc
clear all
close all

clim_Hi  = [0,3500];
clim_phi = [0,7];

%% Read data
foldernames_perfect = {...
  'BIVMIP_A_perfect_40km',...
  'BIVMIP_A_perfect_20km',...
  'BIVMIP_A_perfect_10km'};
foldernames_inverted = {...
  'BIVMIP_A_inv_40km_perfect',...
  'BIVMIP_A_inv_20km_perfect',...
  'BIVMIP_A_inv_10km_perfect'};

for fi = 1: length(foldernames_perfect)
  
  filename_restart       = [foldernames_perfect{fi} '/restart_ANT.nc'];
  filename_help_fields   = [foldernames_perfect{fi} '/help_fields_ANT.nc'];
  
  perfect(fi).x          = ncread(filename_restart    ,'x');
  perfect(fi).y          = ncread(filename_restart    ,'y');
  perfect(fi).time       = ncread(filename_restart    ,'time');
  ti = length(perfect(fi).time);
  
  perfect(fi).Hi         = ncread(filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  perfect(fi).phi_fric   = ncread(filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
end

for fi = 1: length(foldernames_inverted)
  
  filename_restart       = [foldernames_inverted{fi} '/restart_ANT.nc'];
  filename_help_fields   = [foldernames_inverted{fi} '/help_fields_ANT.nc'];
  
  inverted(fi).x         = ncread(filename_restart    ,'x');
  inverted(fi).y         = ncread(filename_restart    ,'y');
  inverted(fi).time      = ncread(filename_restart    ,'time');
  ti = length(inverted(fi).time);
  
  inverted(fi).Hi        = ncread(filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  inverted(fi).phi_fric  = ncread(filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
end

%% Set up GUI
wa = 300;
ha = 300;

margins_hor = [75,25,25,25,120];
margins_ver = [25,50,50];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = zeros( nay,nax);
H.Axa = zeros( nay,nax);
for i = 1: nay
  for j = 1: nax
    x = sum(margins_hor(1:j )) + (j -1)*wa;
    ip = nay+1-i;
    y = sum(margins_ver(1:ip)) + (ip-1)*ha;
    H.Ax( i,j) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],...
      'fontsize',24,'xlim',[-700e3,700e3],'xtick',[],'ylim',[-700e3,700e3],'ytick',[]);
    if (i==1)
      set(H.Ax(i,j),'xaxislocation','top')
    end
  end
end

ylabel(H.Ax(1,1),['Till friction angle (' char(176) ')']);
ylabel(H.Ax(2,1),'Ice thickness (m)');
xlabel(H.Ax(1,1),'Perfect (10 km)');
xlabel(H.Ax(1,2),'Inverted (40 km)');
xlabel(H.Ax(1,3),'Inverted (20 km)');
xlabel(H.Ax(1,4),'Inverted (10 km)');

% Colorbars - right
pos = get(H.Ax(1,4),'position');
H.Cbar14 = colorbar(H.Ax(1,4),'location','eastoutside');
set(H.Ax(1,j),'position',pos);
for j = 1:4
  colormap(H.Ax(1,j),parula(16));
  set(H.Ax(1,j),'clim',clim_phi);
end
ylabel(H.Cbar14,['Till friction angle (' char(176) ')']);

pos = get(H.Ax(2,4),'position');
H.Cbar24 = colorbar(H.Ax(2,4),'location','eastoutside');
set(H.Ax(2,4),'position',pos);
for j = 1:4
  colormap(H.Ax(2,j),itmap(16));
  set(H.Ax(2,j),'clim',clim_Hi);
end
ylabel(H.Cbar24,'Ice thickness (m)');

%% Plot data

% Plot perfect fields
image('parent',H.Ax(1,1),'xdata',perfect(3).x,'ydata',perfect(3).y,'cdata',perfect(3).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(2,1),'xdata',perfect(3).x,'ydata',perfect(3).y,'cdata',perfect(3).Hi','cdatamapping','scaled');

% Plot inverted fields
for j = 1:3
  
  image('parent',H.Ax(1,j+1),'xdata',inverted(j).x,'ydata',inverted(j).y,'cdata',inverted(j).phi_fric','cdatamapping','scaled');
  image('parent',H.Ax(2,j+1),'xdata',inverted(j).x,'ydata',inverted(j).y,'cdata',inverted(j).Hi','cdatamapping','scaled');
  
  % Mark ice-free area as uncertain
  image('parent',H.Ax(1,j+1),'xdata',inverted(j).x,'ydata',inverted(j).y,'cdata',ones(length(inverted(j).x),length(inverted(j).y),3),...
    'alphadata',double(inverted(j).Hi'<0.1)*0.5,'alphadatamapping','none');
end
