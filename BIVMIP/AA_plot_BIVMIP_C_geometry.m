clc
clear all
close all

foldername = '../../BIVMIP/BIVMIP_C_perfect_5km_tuneA';

cmap_phi = parula(256);
clim_phi = [0,6];

filename_restart     = [foldername '/restart_ANT.nc'];
filename_help_fields = [foldername '/help_fields_ANT.nc'];

x        = ncread( filename_restart    ,'x');
y        = ncread( filename_restart    ,'y');
time     = ncread( filename_restart    ,'time'); ti = length(time);
Hi       = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
Hb       = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
Hs       = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
phi_fric = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);

ns = 5;
x        = x(        1:ns:end  );
Hi       = Hi(       1:ns:end,:);
Hb       = Hb(       1:ns:end,:);
Hs       = Hs(       1:ns:end,:);
phi_fric = phi_fric( 1:ns:end,:);

[xg,yg] = meshgrid(x,y);
imid = round(length(x)/2);
jmid = round(length(y)/2);
iend = 1; while Hi(iend,jmid)>0; iend = iend+1; end; iend = iend+1;

Hib = Hs - Hi;
Hib( abs(Hib - Hb) < 1) = NaN;

%% Set up GUI
wa = 700;
ha = 405;

margins_hor = [25,100];
margins_ver = [25,25];

wf = margins_hor(1) + wa + margins_hor(2);
hf = margins_ver(1) + ha + margins_ver(2);

H.Fig = figure('position',[300,300,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'fontsize',24,'xlim',[min(x),max(x)],'ylim',[min(y),max(y)],'zlim',[min(Hb(:)),max(Hs(:))],...
  'xtick',[],'ytick',[],'ztick',[],'cameraposition',[3.3456   -0.4535    0.0186]*1e6);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';
H.Ax.ZAxis.Visible = 'off';

% Bedrock patch with bed roughness colours
zdata = Hb;
surface('parent',H.Ax,'xdata',x,'ydata',y,'zdata',zdata','cdata',phi_fric','edgecolor','k');

% Ice surface patch
zdata = Hs;
zdata(:,1:jmid-1) = NaN;
zdata(iend:end,:) = NaN;
surface('parent',H.Ax,'xdata',x,'ydata',y,'zdata',zdata','facecolor','w','edgecolor','k');

% Ice base patch
zdata = Hib;
zdata(:,1:jmid-1) = NaN;
zdata(iend:end,:) = NaN;
surface('parent',H.Ax,'xdata',x,'ydata',y,'zdata',zdata','facecolor','w','edgecolor','k');

% Nice black line for the surface edge
zdata = Hs(:,jmid);
zdata(iend:end) = NaN;
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',zdata,'color','k','linewidth',3);
zdata = Hs(:,jmid) - Hi(:,jmid);
zdata(iend:end) = NaN;
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',zdata,'color','k','linewidth',3);
zdata = Hb(:,jmid);
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',zdata,'color','k','linewidth',3);

% Colormap for bed roughness
colormap(H.Ax,cmap_phi);
set(H.Ax,'clim',clim_phi);

% Colorbar
pos = get(H.Ax,'position');
pos(4) = pos(4);
H.Ax2 = axes('parent',H.Fig,'units','pixels','position',pos,'fontsize',24,'color','none','clim',clim_phi);
colormap(H.Ax2,cmap_phi);
H.Ax2.XAxis.Visible = 'off';
H.Ax2.YAxis.Visible = 'off';
H.Cbar = colorbar(H.Ax2,'location','eastoutside');
set(H.Ax2,'position',pos);
ylabel(H.Cbar,['Till friction angle (' char(176) ')']);