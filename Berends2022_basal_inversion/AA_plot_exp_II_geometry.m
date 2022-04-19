clc
clear all
close all

foldername = 'exp_II_target_5km';

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

margins_hor = [100,300];
margins_ver = [25,100];

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
H.Cbar = colorbar(H.Ax2,'location','westoutside');
set(H.Ax2,'position',pos);
ylabel(H.Cbar,['Till friction angle (' char(176) ')']);
pos = get(H.Cbar,'position');
pos(1) = 0.07;
set(H.Cbar,'position',pos);

%% Second axis for transect
wa2 = 380;
ha2 = 240;
H.Ax2 = axes('units','pixels','position',[margins_hor(1)+wa/2+250,margins_ver(1)+ha/2+45,wa2,ha2],...
  'fontsize',24,'xlim',[0,800],'ylim',[-1000,3000],'xgrid','on','ygrid','on');
xlabel(H.Ax2,'x (km)');
ylabel(H.Ax2,'z (m)');

x = ncread(filename_restart,'x');
y = ncread(filename_restart,'y');
jmid = find(y==0);
Hi  = ncread(filename_restart,'Hi',[1,jmid,ti],[Inf,1,1]);
Hb  = ncread(filename_restart,'Hb',[1,jmid,ti],[Inf,1,1]);
Hs  = ncread(filename_restart,'Hs',[1,jmid,ti],[Inf,1,1]);
Hib = Hs - Hi;
line('parent',H.Ax2,'xdata',400+x/1e3,'ydata',Hb ,'linewidth',3);
line('parent',H.Ax2,'xdata',400+x/1e3,'ydata',Hs ,'linewidth',3);
line('parent',H.Ax2,'xdata',400+x/1e3,'ydata',Hib,'linewidth',3);