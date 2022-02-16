clc
clear all
close all

%% Read data
filename  = 'BIVMIP_C_perfect_5km/restart_ANT.nc';
filename2 = 'BIVMIP_C_perfect_5km/help_fields_ANT.nc';

time = ncread(filename,'time');
ti = length(time);

Hi = ncread(filename,'Hi',[1,1,ti],[Inf,Inf,1]);
Hb = ncread(filename,'Hb',[1,1,ti],[Inf,Inf,1]);
Hs = Hi + max( -910 / 1028 * Hi, Hb);
Hib = Hs - Hi;
TAF = Hi - max(0, (0 - Hb) * (1028 / 910));

x = ncread(filename,'x'); x = (x + 4e5) / 1e3 - 5;
y = ncread(filename,'y'); y = y / 1e3 - 5;

phi_fric = ncread(filename2,'phi_fric',[1,1,ti],[Inf,Inf,1]);

%% GUI
wa = 800;
ha = 600;

wf = 105 + wa + 55;
hf = 25 + ha + 25;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[-20,-50,wa+25+100,ha+25+50],'fontsize',24);

%% Plot

hold on

% Bedrock
color = [0.6,0.4,0.0];
[ydata,xdata] = meshgrid(y,x);
zdata = Hb;
H.bed = surf('parent',H.Ax,'xdata',xdata,'ydata',ydata,'zdata',zdata,'cdata',phi_fric,'facecolor','interp','edgecolor','none');
for i = 1:4:length(x)
  line('parent',H.Ax,'xdata',zeros(size(y))+x(i),'ydata',y,'zdata',Hb(i,:),'color',color*0.7);
end
for j = 1:length(y)
  line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x))+y(j),'zdata',Hb(:,j),'color',color*0.7);
end
set(H.Ax,'clim',[0,5]);

% Ice surface
color = [0.2,0.4,0.8];
[ydata,xdata] = meshgrid(y,x);
zdata = Hs;
zdata( Hs == min(Hs(:))) = NaN;
i1 = length(x); while isnan(zdata(i1,1)); i1=i1-1; end; zdata(i1,:) = 0;
H.icesurf = surf('parent',H.Ax,'xdata',xdata,'ydata',ydata,'zdata',zdata,'facecolor',color,'facealpha',0.4,'edgecolor','none');
for i = 1:4:i1-1
  zdata = Hs(i,:);
  zdata( zdata == min(Hs(:))) = NaN;
  line('parent',H.Ax,'xdata',zeros(size(y))+x(i),'ydata',y,'zdata',zdata,'color',color*0.7);
end
for j = 1:j
  zdata = Hs(:,j);
  zdata( zdata == min(Hs(:))) = NaN;
  zdata(i1) = 0;
  line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x))+y(j),'zdata',zdata,'color',color*0.7);
end

% Ice base
% color = [0.1,0.8,0.4];
[ydata,xdata] = meshgrid(y,x);
zdata = Hib;
zdata( abs(Hib-Hb)<1.0) = NaN;
zdata( Hs == min(Hs(:))) = NaN;
zdata(i1,:) = 0;
surf('parent',H.Ax,'xdata',xdata,'ydata',ydata,'zdata',zdata,'facecolor',color,'facealpha',0.4,'edgecolor','none');
henk = Hib;
henk( abs(Hib-Hb)<1.0) = NaN;
for i = 1:4:i1-1
  zdata = henk(i,:);
  zdata( zdata == min(Hs(:))) = NaN;
  line('parent',H.Ax,'xdata',zeros(size(y))+x(i),'ydata',y,'zdata',zdata,'color',color*0.7);
end
for j = 1:j
  zdata = henk(:,j);
  zdata( zdata == min(Hs(:))) = NaN;
  zdata(i1) = 0;
  line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x))+y(j),'zdata',zdata,'color',color*0.7);
end

% Make axes prettier
set(H.Ax,'xtick',[],'ytick',[],'ztick',[],'xlim',[0,800],'ylim',[-120,120],'zlim',[-800,3000],...
  'cameraposition',[0.4153   -0.1214    2.0985]*1e4);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';
H.Ax.ZAxis.Visible = 'off';

% Grounding line
x_GL = zeros(size(y));
for j = 1:length(y)
  i = 1;
  while TAF(i)>0; i=i+1; end
  TAF1 = TAF(i-1,j);
  TAF2 = TAF(i,j);
  lambda = TAF1 / (TAF1 - TAF2);
  x_GL(j) = (1 - lambda) * x(i-1) + lambda * x(i);
end

xdata = x_GL;
ydata = y;

% On bedrock
zdata = zeros(size(y));
for j = 1:length(y)
  yq = y(j);
  xq = x_GL(j);
  zdata(j) = interp2(x,y,Hb',xq,yq);
end
line('parent',H.Ax,'xdata',xdata,'ydata',ydata,'zdata',zdata,'linewidth',5,'color','r');

% On ice surface
zdata = zeros(size(y));
for j = 1:length(y)
  yq = y(j);
  xq = x_GL(j);
  zdata(j) = interp2(x,y,Hs',xq,yq);
end
line('parent',H.Ax,'xdata',xdata,'ydata',ydata,'zdata',zdata,'linewidth',5,'color','r');

%% Colorbar
pos = get(H.Ax,'position');
H.Cbar = colorbar(H.Ax,'location','westoutside');
set(H.Ax,'position',pos);
pos = get(H.Cbar,'position');
pos(1) = 0.08;
pos(2) = 0.1;
pos(4) = 0.8;
set(H.Cbar,'position',pos);
ylabel(H.Cbar,['Till friction angle (' char(176) ')']);

%% Transect
% Line in 3D model
jmid = find(y==0);
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',Hb(:,jmid),'linewidth',5,'color','k');
zdata = Hs;
zdata(Hi==0) = NaN;
zdata(i1,:) = 0;
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',zdata(:,jmid),'linewidth',5,'color','k');
zdata = Hib;
zdata(Hi==0) = NaN;
zdata(i1,:) = 0;
line('parent',H.Ax,'xdata',x,'ydata',zeros(size(x)),'zdata',zdata(:,jmid),'linewidth',5,'color','k');

% Actual cross-section
fx = 0.5;
fy = 0.4;
H.Ax2 = axes('parent',H.Fig,'units','pixels','position',[125+wa*(1-fx),25+ha*(1-fy),wa*fx,ha*fy],...
  'xlim',[0,800],'ylim',[-800,3000],'fontsize',24,'xgrid','on','ygrid','on');
xlabel(H.Ax2,'x (km)')
ylabel(H.Ax2,'z (m)')
jmid = ceil(length(y)/2);
line('parent',H.Ax2,'xdata',x,'ydata',Hb(:,jmid),'linewidth',2,'color','k');
zdata = Hs;
zdata(Hi==0) = NaN;
zdata(i1,:) = 0;
line('parent',H.Ax2,'xdata',x,'ydata',zdata(:,jmid),'linewidth',2,'color','k');
zdata = Hib;
zdata(Hi==0) = NaN;
zdata(i1,:) = 0;
line('parent',H.Ax2,'xdata',x,'ydata',zdata(:,jmid),'linewidth',2,'color','k');