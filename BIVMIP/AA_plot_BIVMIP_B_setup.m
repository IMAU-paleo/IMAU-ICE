clc
clear all
close all

filename_restart     = 'BIVMIP_B_perfect_10km/restart_ANT.nc';
filename_help_fields = 'BIVMIP_B_perfect_10km/help_fields_ANT.nc';

% Read data
x        = ncread( filename_restart    ,'x'); nx = length(x);
y        = ncread( filename_restart    ,'y'); ny = length(y);
time     = ncread( filename_restart    ,'time');
ti = length(time);
Hi       = ncread( filename_restart    ,'Hi'      ,[1,1,  ti],[Inf,Inf,  1]);
Hb       = ncread( filename_restart    ,'Hb'      ,[1,1,  ti],[Inf,Inf,  1]);
Hs       = ncread( filename_restart    ,'Hs'      ,[1,1,  ti],[Inf,Inf,  1]);
phi_fric = ncread( filename_help_fields,'phi_fric',[1,1,  ti],[Inf,Inf,  1]);
u_surf   = ncread( filename_help_fields,'u_3D'    ,[1,1,1,ti],[Inf,Inf,1,1]);
v_surf   = ncread( filename_help_fields,'v_3D'    ,[1,1,1,ti],[Inf,Inf,1,1]);

% Plot stuff
wa = 500;
ha = 350;

wf = 25 + wa + 80;
hf = 25 + ha + 25;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax1  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],...
  'xlim',[min(x),max(x)],'ylim',[min(y),max(y)],'fontsize',24,'color',[0.5,0.5,1.0],...
  'xtick',[],'ytick',[],'ztick',[]);

% Till friction angle as bedrock
image('parent',H.Ax1,'xdata',x,'ydata',y,'cdata',phi_fric,'cdatamapping','scaled')

% Ice sheet in 3D
[xg,yg] = meshgrid(x,y);
z = Hs - 0.1;
z( yg < 0) = NaN;
surface('parent',H.Ax1,'xdata',xg,'ydata',yg,'zdata',z,'facecolor',0.99*[1,1,1],'edgecolor','none','facealpha',1,...
  'facelighting','gouraud','ambientstrength',0.5,'diffusestrength',0.9,'specularstrength',0.3);
light('parent',H.Ax1,'style','local','position',[max(x)/2,min(y)*3/4,5000]);
set(H.Ax1,'cameraposition',1e6 * [-6.2035   -4.8662    0.0225]);

% Mesh for better shape recognition
i0 = find(x==0);
Hsi0 = Hs(i0,:);
line('parent',H.Ax1,'ydata',zeros(size(y))+y(i0),'xdata',x,'zdata',Hsi0,'color','k','linewidth',2);
for i = i0:6:length(x)
  Hsi = Hs(i,:);
  line('parent',H.Ax1,'ydata',zeros(size(y))+x(i),'xdata',y,'zdata',Hsi,'color','k');
end
for j = 1:6:length(y)
  Hsi = Hs(:,j);
  Hsi(y<0) = NaN;
  line('parent',H.Ax1,'xdata',zeros(size(y))+y(j),'ydata',y,'zdata',Hsi,'color','k');
end

% Colorbar
H.Cbar1 = colorbar(H.Ax1,'location','eastoutside');
ylabel(H.Cbar1,['\phi_{fric} (' char(176) ')']);