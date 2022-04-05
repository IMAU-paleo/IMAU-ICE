clc
clear all
close all

foldernames = {...
  'exp_I_target_10km',...
  'exp_I_inv_10km_unperturbed',...
  'exp_I_inv_10km_visc_hi',...
  'exp_I_inv_10km_visc_lo',...
  'exp_I_inv_10km_SMB_hi',...
  'exp_I_inv_10km_SMB_lo',...
  'exp_I_inv_10km_p_hi',...
  'exp_I_inv_10km_p_lo',...
  'exp_I_inv_10km_ut_hi',...
  'exp_I_inv_10km_ut_lo'};

cmap_phi = parula(16);
cmap_Hs  = itmap(16);

clim_phi = [0,10];
clim_Hs  = [0,3000];

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

wa = 200;
ha = 200;

margins_hor = [25,25,25,25,25,25,130];
margins_ver = [25,25,25,25,50];

margins_hor = [25,5,5,5,5,5,130];
margins_ver = [25,5,5,5,50];

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

set(H.Ax(2,1),'visible','off')
set(H.Ax(2,2),'visible','off')
set(H.Ax(4,1),'visible','off')
set(H.Ax(4,2),'visible','off')

for j = 1:6
  colormap(H.Ax(1,j),cmap_phi);
  colormap(H.Ax(2,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi);
  set(H.Ax(2,j),'clim',clim_phi);

  colormap(H.Ax(3,j),cmap_Hs);
  colormap(H.Ax(4,j),cmap_Hs);
  set(H.Ax(3,j),'clim',clim_Hs);
  set(H.Ax(4,j),'clim',clim_Hs);
end

xlabel( H.Ax(1,1),'Target');
xlabel( H.Ax(1,2),'Unperturbed');
xlabel( H.Ax(1,3),'Viscosity');
xlabel( H.Ax(1,4),'SMB');
xlabel( H.Ax(1,5),'Z-I p');
xlabel( H.Ax(1,6),'Z-I u\_t');

% Colorbars
pos1 = get(H.Ax(1,6),'position');
pos2 = get(H.Ax(2,6),'position');
xlo = pos1(1)+pos1(3)+25;
xhi = xlo + 125;
ylo = pos2(2);
yhi = pos1(2)+pos1(4);
H.Axcbar1 = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
H.Axcbar1.XAxis.Visible = 'off';
H.Axcbar1.YAxis.Visible = 'off';
set(H.Axcbar1,'clim',clim_phi);
colormap(H.Axcbar1,cmap_phi);
pos = get(H.Axcbar1,'position');
H.Cbar1 = colorbar(H.Axcbar1,'location','west');
set(H.Axcbar1,'position',pos);
ylabel(H.Cbar1,['Till friction angle (' char(176) ')']);

pos1 = get(H.Ax(3,6),'position');
pos2 = get(H.Ax(4,6),'position');
xlo = pos1(1)+pos1(3)+25;
xhi = xlo + 125;
ylo = pos2(2);
yhi = pos1(2)+pos1(4);
H.Axcbar2 = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
H.Axcbar2.XAxis.Visible = 'off';
H.Axcbar2.YAxis.Visible = 'off';
set(H.Axcbar2,'clim',clim_Hs);
colormap(H.Axcbar2,cmap_Hs);
pos = get(H.Axcbar2,'position');
H.Cbar2 = colorbar(H.Axcbar2,'location','west');
set(H.Axcbar2,'position',pos);
ylabel(H.Cbar2,'Surface elevation (m)');

%% Plot results

for i = 1:4
  for j = 1:6
    x = results(j).x;
    y = results(j).y;
    set(H.Ax(i,j),'xlim',[min(x),max(x)]*0.9,'ylim',[min(y),max(y)]*0.9);

    % Blank white image to cover the axes lines (because they're ugly)
    cdata = zeros(length(y),length(x),3)+1;
    image('parent',H.Ax(i,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
    image('parent',H.Ax(i,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  end
end
  
% Top two rows: till friction angle
ax = H.Ax(1,1);
R  = results(1);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(1,2);
R  = results(2);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(1,3);
R  = results(3);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(1,4);
R  = results(5);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(1,5);
R  = results(7);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(1,6);
R  = results(9);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,3);
R  = results(4);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,4);
R  = results(6);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,5);
R  = results(8);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,6);
R  = results(10);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Bottom two rows: surface elevation
ax = H.Ax(3,1);
R  = results(1);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(3,2);
R  = results(2);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(3,3);
R  = results(3);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(3,4);
R  = results(5);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(3,5);
R  = results(7);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(3,6);
R  = results(9);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,3);
R  = results(4);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,4);
R  = results(6);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,5);
R  = results(8);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,6);
R  = results(10);
cdata = R.Hs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);