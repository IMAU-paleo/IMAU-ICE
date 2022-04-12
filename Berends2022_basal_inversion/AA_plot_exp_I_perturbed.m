clc
clear all
close all

foldername_target = 'exp_I_target_20km';

foldernames = {...
  'exp_I_inv_20km_visc_hi',...          %  1
  'exp_I_inv_20km_visc_lo',...          %  2
  'exp_I_inv_20km_SMB_hi',...           %  3
  'exp_I_inv_20km_SMB_lo',...           %  4
  'exp_I_inv_20km_topo_hi',...          %  5
  'exp_I_inv_20km_topo_lo',...          %  6
  'exp_I_inv_20km_p_hi',...             %  7
  'exp_I_inv_20km_p_lo',...             %  8
  'exp_I_inv_20km_ut_hi',...            %  9
  'exp_I_inv_20km_ut_lo'};              % 10

cmap_phi = jet(32);
cmap_Hs  = flipud(lbmap(32,'redblue'));

clim_phi = [-8,8];
clim_Hs  = [-120,120];

%% Read data

% Target
filename_restart     = [foldername_target '/restart_ANT.nc'];
filename_help_fields = [foldername_target '/help_fields_ANT.nc'];

target.x         = ncread( filename_restart    ,'x');
target.y         = ncread( filename_restart    ,'y');
target.time      = ncread( filename_restart    ,'time'); ti = length(target.time);
target.Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
target.Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
target.Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
target.phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);

% Perturbed runs
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
  
  results(fi).dphi_fric = results(fi).phi_fric - target.phi_fric;
  results(fi).dHs       = results(fi).Hs       - target.Hs;
end

%% Set up GUI

wa = 180;
ha = 180;

margins_hor = [45,5,5,5,5,130];
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

for j = 1:nax
  colormap(H.Ax(1,j),cmap_phi);
  colormap(H.Ax(2,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi);
  set(H.Ax(2,j),'clim',clim_phi);

  colormap(H.Ax(3,j),cmap_Hs);
  colormap(H.Ax(4,j),cmap_Hs);
  set(H.Ax(3,j),'clim',clim_Hs);
  set(H.Ax(4,j),'clim',clim_Hs);
end

xlabel( H.Ax(1,1),'Viscosity');
xlabel( H.Ax(1,2),'SMB');
xlabel( H.Ax(1,3),'Topography');
xlabel( H.Ax(1,4),'Z-I p');
xlabel( H.Ax(1,5),'Z-I u\_t');

ylabel( H.Ax(1,1),'High')
ylabel( H.Ax(2,1),'Low')
ylabel( H.Ax(3,1),'High')
ylabel( H.Ax(4,1),'Low')

% Colorbars
pos1 = get(H.Ax(1,nax),'position');
pos2 = get(H.Ax(2,nax),'position');
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
ylabel(H.Cbar1,['\Delta till friction angle (' char(176) ')']);

pos1 = get(H.Ax(3,nax),'position');
pos2 = get(H.Ax(4,nax),'position');
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
ylabel(H.Cbar2,'\Delta surface elevation (m)');

%% Plot results

for i = 1:4
  for j = 1:nax
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

% Visc
ax = H.Ax(1,1);
R  = results(1);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,1);
R  = results(2);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% SMB
ax = H.Ax(1,2);
R  = results(3);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,2);
R  = results(4);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Topo
ax = H.Ax(1,3);
R  = results(5);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,3);
R  = results(6);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% ut
ax = H.Ax(1,4);
R  = results(7);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,4);
R  = results(8);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% p
ax = H.Ax(1,5);
R  = results(9);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(2,5);
R  = results(10);
cdata = R.dphi_fric';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Bottom two rows: surface elevation

% Viscosity
ax = H.Ax(3,1);
R  = results(1);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,1);
R  = results(2);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% SMB
ax = H.Ax(3,2);
R  = results(3);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,2);
R  = results(4);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Topo
ax = H.Ax(3,3);
R  = results(5);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,3);
R  = results(6);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% ut
ax = H.Ax(3,4);
R  = results(7);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,4);
R  = results(8);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% p
ax = H.Ax(3,5);
R  = results(9);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

ax = H.Ax(4,5);
R  = results(10);
cdata = R.dHs';
adata = zeros(size(cdata));
adata( R.Hi'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);