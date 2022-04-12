clc
clear all
close all

foldernames = {...
  'exp_II_target_5km',...
  'exp_II_target_5km',...
  'exp_II_inv_5km_unperturbed',...
  'exp_II_inv_5km_unperturbed'};

cmap_phi  = parula(256);
cmap_Hs   = itmap(16);

cmap_dphi = jet(32);
cmap_dHs  = flipud(lbmap(32,'redblue'));

clim_phi  = [0,6];
clim_Hs   = [0,2700];

clim_dphi = [-8,8];
clim_dHs  = [-120,120];

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
  
  % Remove the weird artefact in the northwest corner
  results(fi).Hi(1,:) = min(results(fi).Hi(1,:));
  results(fi).Hs(1,:) = min(results(fi).Hs(1,:));
  
  % Calculate thickness above flotation
  ice_density           =  910.0;
  seawater_density      = 1028.0;
  results(fi).TAF       = results(fi).Hi - max(0, (-results(fi).Hb) * (seawater_density / ice_density));
  
end

% Calculate errors
for fi = 3:4
  results(fi).dHs       = results(fi).Hs       - results(fi-2).Hs;
  results(fi).dphi_fric = results(fi).phi_fric - results(fi-2).phi_fric;
end

%% Set up GUI

wa = 400;
ha = 100;

margins_hor = [25,25,25];
margins_ver = [90,50,120,50];

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

colormap(H.Ax(1,1),cmap_phi);
set(H.Ax(1,1),'clim',clim_phi);

colormap(H.Ax(1,2),cmap_Hs);
set(H.Ax(1,2),'clim',clim_Hs);

for j = 2:nay
  colormap(H.Ax(j,1),cmap_dphi);
  set(H.Ax(j,1),'clim',clim_dphi);
  
  colormap(H.Ax(j,2),cmap_dHs);
  set(H.Ax(j,2),'clim',clim_dHs);
end

lab = xlabel(H.Ax(1,1),'Target');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

lab = xlabel(H.Ax(2,1),'Unperturbed (5 km)');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

lab = xlabel(H.Ax(3,1),'Unperturbed (2 km)');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

x = results(2).x;
y = results(2).y;
for i = 1:nay
  for j = 1:nax
    set(H.Ax(i,j),'xlim',[min(x),max(x)*0.5],'ylim',[min(y),max(y)]);
  end
end

%% Colorbars
xlo = margins_hor(1);
xhi = xlo+wa;
ylo = margins_ver(1) + ha + margins_ver(2) + ha + 5;
yhi = margins_ver(1) + ha + margins_ver(2) + ha + margins_ver(3) - 5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24,'color','none');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_phi);
set(ax,'clim',clim_phi);
H.Cbar1 = colorbar(ax,'location','north');
xlabel(H.Cbar1,['Till friction angle (' char(176) ')']);

xlo = margins_hor(1)+wa+margins_hor(2);
xhi = xlo+wa;
ylo = margins_ver(1) + ha + margins_ver(2) + ha + 5;
yhi = margins_ver(1) + ha + margins_ver(2) + ha + margins_ver(3) - 5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24,'color','none');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_Hs);
set(ax,'clim',clim_Hs);
H.Cbar2 = colorbar(ax,'location','north');
xlabel(H.Cbar2,'Surface elevation (m)');

xlo = margins_hor(1);
xhi = xlo+wa;
ylo = 5;
yhi = margins_ver(1)-5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_dphi);
set(ax,'clim',clim_dphi);
H.Cbar1 = colorbar(ax,'location','north');
xlabel(H.Cbar1,['\Delta till friction angle (' char(176) ')']);

xlo = margins_hor(1)+wa+margins_hor(2);
xhi = xlo+wa;
ylo = 5;
yhi = margins_ver(1)-5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_dHs);
set(ax,'clim',clim_dHs);
H.Cbar2 = colorbar(ax,'location','north');
xlabel(H.Cbar2,'\Delta surface elevation (m)');

%% Plot results - target
  
% Construct target grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL_target = contour('parent',H.tempax,'xdata',results(1).y,'ydata',results(1).x,'zdata',results(1).TAF,'levellist',0);
close(H.tempfig);

% Blank white image to cover the axes lines (because they're ugly)
cdata = zeros(2,2,3)+1;
image('parent',H.Ax(1,1),'xdata',get(H.Ax(1,1),'xlim')*1.1,'ydata',get(H.Ax(1,1),'ylim')*1.1,'cdata',cdata);
image('parent',H.Ax(1,2),'xdata',get(H.Ax(1,2),'xlim')*1.1,'ydata',get(H.Ax(1,1),'ylim')*1.1,'cdata',cdata);

% Left column: till friction angle
ax = H.Ax(1,1);
R  = results(2);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.TAF'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Right column: surface elevation
ax = H.Ax(1,2);
R  = results(2);
cdata = R.Hs';
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

% Add target grounding line contour
C = C_GL_target;
while ~isempty(C)
  n  = C(2,1);
  Ct = C(:,2:2+n-1);
  C = C(:,2+n:end);
  line('parent',H.Ax(1,1),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
  line('parent',H.Ax(1,2),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
end

%% Plot results - inverted

for j = 1:2
  
  jr = j+2;
  ja = j+1;

  % Left column: till friction angle
  ax = H.Ax(ja,1);
  R  = results(jr);
  cdata = R.dphi_fric';
  adata = zeros(size(cdata));
  adata( R.TAF'>0) = 1;
  image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

  % Right column: surface elevation
  ax = H.Ax(ja,2);
  R  = results(jr);
  cdata = R.dHs';
  image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

  % Add target grounding line contour
  C = C_GL_target;
  while ~isempty(C)
    n  = C(2,1);
    Ct = C(:,2:2+n-1);
    C = C(:,2+n:end);
    line('parent',H.Ax(ja,1),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
    line('parent',H.Ax(ja,2),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
  end
  
  % Construct inverted-geometry grounding line contour
  H.tempfig = figure;
  H.tempax  = axes('parent',H.tempfig);
  C_GL = contour('parent',H.tempax,'xdata',results(jr).y,'ydata',results(jr).x,'zdata',results(jr).TAF,'levellist',0);
  close(H.tempfig);

  % Add inverted-geometry grounding line contour
  C = C_GL;
  while ~isempty(C)
    n  = C(2,1);
    Ct = C(:,2:2+n-1);
    C = C(:,2+n:end);
    line('parent',H.Ax(ja,2),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',2,'color','k','linestyle','--');
    line('parent',H.Ax(ja,1),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',2,'color','k','linestyle','--');
  end
  
end