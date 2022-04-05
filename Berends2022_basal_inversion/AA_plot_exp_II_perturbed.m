clc
clear all
close all

foldernames = {...
  'BIVMIP_C_retreat_2km_perfect',...
  'BIVMIP_C_retreat_2km_visc_hi',...
  'BIVMIP_C_retreat_2km_visc_lo',...
  'BIVMIP_C_retreat_2km_SMB_hi',...
  'BIVMIP_C_retreat_2km_SMB_lo'};

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
  time      = ncread( filename_restart    ,'time'); ti = 1;
  Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
  Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
  phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  
  ice_density                      =  910.0;
  seawater_density                 = 1028.0;
  TAF       = Hi - max(0, (-Hb) * (seawater_density / ice_density));
  
  results(fi).x        = x;
  results(fi).y        = y;
  results(fi).Hi       = Hi;
  results(fi).Hb       = Hb;
  results(fi).Hs       = Hs;
  results(fi).phi_fric = phi_fric;
  results(fi).TAF      = TAF;
end

%% Set up GUI

wa = 400;
ha = 100;

margins_hor = [50,25,25];
margins_ver = [90,5,30,5,30,15];

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

for j = 1:nay
  colormap(H.Ax(j,1),cmap_phi);
  set(H.Ax(j,1),'clim',clim_phi);
  
  colormap(H.Ax(j,2),cmap_Hs);
  set(H.Ax(j,2),'clim',clim_Hs);
end

%% "invisible" axeses for y-labels
pos1 = get(H.Ax(1,1),'position');
pos2 = get(H.Ax(1,1),'position');
xlo = pos1(1);
xhi = xlo;
ylo = pos2(2);
yhi = pos1(2)+pos1(4);
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24,'xlim',[0,1],'ylim',[0,1],...
  'xtick',[],'ytick',[]);
cdata = zeros(2,2,3)+1;
image('parent',ax,'xdata',[-1,2],'ydata',[-1,2],'cdata',cdata);
ylabel(ax,'Perfect');

pos1 = get(H.Ax(2,1),'position');
pos2 = get(H.Ax(3,1),'position');
xlo = pos1(1);
xhi = xlo;
ylo = pos2(2);
yhi = pos1(2)+pos1(4);
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24,'xlim',[0,1],'ylim',[0,1],...
  'xtick',[],'ytick',[]);
cdata = zeros(2,2,3)+1;
image('parent',ax,'xdata',[-1,2],'ydata',[-1,2],'cdata',cdata);
ylabel(ax,'Viscosity');

pos1 = get(H.Ax(4,1),'position');
pos2 = get(H.Ax(5,1),'position');
xlo = pos1(1);
xhi = xlo;
ylo = pos2(2);
yhi = pos1(2)+pos1(4);
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24,'xlim',[0,1],'ylim',[0,1],...
  'xtick',[],'ytick',[]);
cdata = zeros(2,2,3)+1;
image('parent',ax,'xdata',[-1,2],'ydata',[-1,2],'cdata',cdata);
ylabel(ax,'SMB');

%% Colorbars
xlo = margins_hor(1);
xhi = xlo+wa;
ylo = 5;
yhi = margins_ver(1)-5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_phi);
set(ax,'clim',clim_phi);
H.Cbar1 = colorbar(ax,'location','north');
xlabel(H.Cbar1,['Till friction angle (' char(176) ')']);

xlo = margins_hor(1)+wa+margins_hor(2);
xhi = xlo+wa;
ylo = 5;
yhi = margins_ver(1)-5;
ax = axes('parent',H.Fig,'units','pixels','position',[xlo,ylo,xhi-xlo,yhi-ylo],'fontsize',24);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
colormap(ax,cmap_Hs);
set(ax,'clim',clim_Hs);
H.Cbar2 = colorbar(ax,'location','north');
xlabel(H.Cbar2,'Surface elevation (m)');

%% Plot results
  
% Construct "perfect" grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL_perfect = contour('parent',H.tempax,'xdata',results(1).y,'ydata',results(1).x,'zdata',results(1).TAF,'levellist',0);
close(H.tempfig);

for i = 1:5
  
  % Construct actual grounding line contour
  H.tempfig = figure;
  H.tempax  = axes('parent',H.tempfig);
  C_GL = contour('parent',H.tempax,'xdata',results(i).y,'ydata',results(i).x,'zdata',results(i).TAF,'levellist',0);
  close(H.tempfig);

  for j = 1:2
    x = results(j).x;
    y = results(j).y;
    set(H.Ax(i,j),'xlim',[min(x),max(x)*0.5],'ylim',[min(y),max(y)]);

    % Blank white image to cover the axes lines (because they're ugly)
    cdata = zeros(length(y),length(x),3)+1;
    image('parent',H.Ax(i,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
    image('parent',H.Ax(i,j),'xdata',x*1.1,'ydata',y*1.1,'cdata',cdata);
    
    % Left column: till friction angle
    if j==1
      ax = H.Ax(i,1);
      R  = results(i);
      cdata = R.phi_fric';
      adata = zeros(size(cdata));
      adata( R.TAF'>0) = 1;
      image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
    end
  
    % Right column: surface elevation
    if j==2
      ax = H.Ax(i,2);
      R  = results(i);
      cdata = R.Hs';
      image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');
    end

    % Add "perfect" grounding line contour
    C = C_GL_perfect;
    while ~isempty(C)
      n  = C(2,1);
      Ct = C(:,2:2+n-1);
      C = C(:,2+n:end);
      line('parent',H.Ax(i,j),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',3,'color','r');
    end

    % Add actual grounding line contour
    C = C_GL;
    while ~isempty(C)
      n  = C(2,1);
      Ct = C(:,2:2+n-1);
      C = C(:,2+n:end);
      line('parent',H.Ax(i,j),'xdata',Ct(2,:),'ydata',Ct(1,:),'linewidth',2,'color','k','linestyle','--');
    end

  end
end
  