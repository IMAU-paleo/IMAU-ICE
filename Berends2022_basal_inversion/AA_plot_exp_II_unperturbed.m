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

cmap_u    = sentinelmap( 32);
cmap_du   = jet( 32);

clim_phi  = [0,6];
clim_Hs   = [0,2700];

clim_dphi = [0.1,10];
clim_dHs  = [-250,250];

clim_u    = [1,1000];
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
  results(fi).dphi_fric = results(fi).phi_fric ./ results(fi-2).phi_fric;
  results(fi).dHs       = results(fi).Hs       -  results(fi-2).Hs;
  results(fi).du        = results(fi).uabs     -  results(fi-2).uabs;
end

%% Set up GUI

wa = 300;
ha = 75;

margins_hor = [25,50,25,25];
margins_ver = [90,90,90,50];

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

for j = 1:1
  colormap(H.Ax(1,j),cmap_phi);
  set(H.Ax(1,j),'clim',clim_phi);
  
  colormap(H.Ax(2,j),cmap_Hs);
  set(H.Ax(2,j),'clim',clim_Hs);
  
  colormap(H.Ax(3,j),cmap_u);
  set(H.Ax(3,j),'clim',clim_u);
end

for j = 2:3
  colormap(H.Ax(1,j),cmap_dphi);
  set(H.Ax(1,j),'clim',clim_dphi,'colorscale','log');
  
  colormap(H.Ax(2,j),cmap_dHs);
  set(H.Ax(2,j),'clim',clim_dHs);
  
  colormap(H.Ax(3,j),cmap_du);
  set(H.Ax(3,j),'clim',clim_du);
end

xlabel( H.Ax(1,1),'Target (5 km)');
xlabel( H.Ax(1,2),'Unperturbed (5 km)');
xlabel( H.Ax(1,3),'Unperturbed (2.5 km)');

x = results(2).x;
y = results(2).y;
for i = 1:nay
  for j = 1:nax
    set(H.Ax(i,j),'xlim',[min(x),max(x)*0.5],'ylim',[min(y),max(y)]);
  end
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
H.Cbarl1 = colorbar(H.Ax(1,1),'location','southoutside');
set(H.Ax(1,1),'position',pos);
ylabel( H.Cbarl1,['Till friction angle (' char(176) ')'])

pos = get(H.Ax(2,1),'position');
H.Cbarl2 = colorbar(H.Ax(2,1),'location','southoutside');
set(H.Ax(2,1),'position',pos);
ylabel( H.Cbarl2,'Surface elevation (m)')

pos = get(H.Ax(3,1),'position');
H.Cbarl3 = colorbar(H.Ax(3,1),'location','southoutside');
set(H.Ax(3,1),'position',pos);
ylabel( H.Cbarl3,'Surface velocity (m yr^{-1})')
set(H.Cbarl3,'ticks',[1,10,100,1000]);

pos1 = get(H.Ax(1,2),'position');
pos2 = get(H.Ax(1,3),'position');
xl = pos1(1);
xu = pos2(1)+pos2(3);
yl = pos1(2);
yu = pos1(2)+pos1(4);
xc = (xl+xu)/2;
xl = xc - wa/2;
xu = xc + wa/2;
ax = axes('parent',H.Fig,'units','pixels','position',[xl,yl,xu-xl,yu-yl],'clim',clim_dphi,'colorscale','log',...
  'fontsize',24,'xtick',[],'ytick',[],'color','none');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_dphi);
pos = get(ax,'position');
H.Cbarr1 = colorbar(ax,'location','southoutside');
set(ax,'position',pos);
ylabel( H.Cbarr1,['\Delta till friction angle (' char(176) ')'])

pos1 = get(H.Ax(2,2),'position');
pos2 = get(H.Ax(2,3),'position');
xl = pos1(1);
xu = pos2(1)+pos2(3);
yl = pos1(2);
yu = pos1(2)+pos1(4);
xc = (xl+xu)/2;
xl = xc - wa/2;
xu = xc + wa/2;
ax = axes('parent',H.Fig,'units','pixels','position',[xl,yl,xu-xl,yu-yl],'clim',clim_dHs,...
  'fontsize',24,'xtick',[],'ytick',[],'color','none');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_dHs);
pos = get(ax,'position');
H.Cbarr2 = colorbar(ax,'location','southoutside');
set(ax,'position',pos);
ylabel( H.Cbarr2,'\Delta surface elevation (m)')

pos1 = get(H.Ax(3,2),'position');
pos2 = get(H.Ax(3,3),'position');
xl = pos1(1);
xu = pos2(1)+pos2(3);
yl = pos1(2);
yu = pos1(2)+pos1(4);
xc = (xl+xu)/2;
xl = xc - wa/2;
xu = xc + wa/2;
ax = axes('parent',H.Fig,'units','pixels','position',[xl,yl,xu-xl,yu-yl],'clim',clim_du,...
  'fontsize',24,'xtick',[],'ytick',[],'color','none');
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_du);
pos = get(ax,'position');
H.Cbarr3 = colorbar(ax,'location','southoutside');
set(ax,'position',pos);
ylabel( H.Cbarr3,'\Delta surface velocity (m yr^{-1})')

%% Plot results - target
  
% Construct target grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL_target = contour('parent',H.tempax,'xdata',results(1).y,'ydata',results(1).x,'zdata',results(1).TAF,'levellist',0);
close(H.tempfig);

% Blank white image to cover the axes lines (because they're ugly)
cdata = zeros(2,2,3)+1;
image('parent',H.Ax(1,1),'xdata',get(H.Ax(1,1),'xlim')*1.1,'ydata',get(H.Ax(1,1),'ylim')*1.1,'cdata',cdata);
image('parent',H.Ax(2,1),'xdata',get(H.Ax(2,1),'xlim')*1.1,'ydata',get(H.Ax(2,1),'ylim')*1.1,'cdata',cdata);
image('parent',H.Ax(3,1),'xdata',get(H.Ax(3,1),'xlim')*1.1,'ydata',get(H.Ax(3,1),'ylim')*1.1,'cdata',cdata);

% Top row: till friction angle
ax = H.Ax(1,1);
R  = results(2);
cdata = R.phi_fric';
adata = zeros(size(cdata));
adata( R.TAF'>0) = 1;
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

% Middle row: surface elevation
ax = H.Ax(2,1);
R  = results(2);
cdata = R.Hs';
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

% Bottom row: surface velocity
ax = H.Ax(3,1);
R  = results(2);
cdata = R.u';
set(ax,'colorscale','log');
image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

% Add target grounding line contour
C = C_GL_target;
while ~isempty(C)
  n  = C(2,1);
  Ct = C(:,2:2+n-1);
  C = C(:,2+n:end);
  line('parent',H.Ax(2,1),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
end

%% Plot results - inverted

for j = 1:2
  
  jr = j+2;
  ja = j+1;
  
  x = results(jr).x;
  y = results(jr).y;

  % Blank white image to cover the axes lines (because they're ugly)
  cdata = zeros(length(y),length(x),3)+1;
  image('parent',H.Ax(1,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  image('parent',H.Ax(2,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
  image('parent',H.Ax(3,ja),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);

  % Top row: till friction angle
  ax = H.Ax(1,ja);
  R  = results(jr);
  cdata = R.dphi_fric';
  adata = zeros(size(cdata));
  adata( R.TAF'>0) = 1;
  image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);

  % Middle row: surface elevation
  ax = H.Ax(2,ja);
  R  = results(jr);
  cdata = R.dHs';
  image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

  % Bottom row: surface velocity
  ax = H.Ax(3,ja);
  R  = results(jr);
  cdata = R.du';
  image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');

  % Add target grounding line contour
  C = C_GL_target;
  while ~isempty(C)
    n  = C(2,1);
    Ct = C(:,2:2+n-1);
    C = C(:,2+n:end);
    line('parent',H.Ax(2,ja),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
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
    line('parent',H.Ax(2,ja),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','k','linestyle','--');
  end
  
end