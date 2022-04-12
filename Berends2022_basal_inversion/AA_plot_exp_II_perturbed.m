clc
clear all
close all

foldername_target = 'exp_II_target_5km';

foldernames = {...
  'exp_II_inv_5km_visc_hi',...          %  1
  'exp_II_inv_5km_visc_lo',...          %  2
  'exp_II_inv_5km_SMB_hi',...           %  3
  'exp_II_inv_5km_SMB_lo',...           %  4
  'exp_II_inv_5km_topo_hi',...          %  5
  'exp_II_inv_5km_topo_lo'};            %  6

cmap_phi = jet(32);
cmap_Hs  = flipud(lbmap(32,'redblue'));

clim_phi = [-8,8];
clim_Hs  = [-120,120];

%% Read data
  
filename_restart     = [foldername_target '/restart_ANT.nc'];
filename_help_fields = [foldername_target '/help_fields_ANT.nc'];

target.x         = ncread( filename_restart    ,'x');
target.y         = ncread( filename_restart    ,'y');
target.time      = ncread( filename_restart    ,'time'); ti = length(target.time);
target.Hi        = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
target.Hb        = ncread( filename_restart    ,'Hb'      ,[1,1,ti],[Inf,Inf,1]);
target.Hs        = ncread( filename_restart    ,'Hs'      ,[1,1,ti],[Inf,Inf,1]);
target.phi_fric  = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  
% Remove the weird artefact in the northwest corner
target.Hi(1,:) = min(target.Hi(1,:));
target.Hs(1,:) = min(target.Hs(1,:));

ice_density      =  910.0;
seawater_density = 1028.0;
target.TAF       = target.Hi - max(0, (-target.Hb) * (seawater_density / ice_density));
  
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
  
  results(fi).dphi_fric = results(fi).phi_fric - target.phi_fric;
  results(fi).dHs       = results(fi).Hs       - target.Hs;
  
  ice_density                      =  910.0;
  seawater_density                 = 1028.0;
  results(fi).TAF       = results(fi).Hi - max(0, (-results(fi).Hb) * (seawater_density / ice_density));
end

%% Set up GUI

wa = 400;
ha = 100;

margins_hor = [50,25,25];
margins_ver = [90,5,50,5,50,5,50];

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

lab = xlabel(H.Ax(1,1),'Viscosity');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

lab = xlabel(H.Ax(3,1),'SMB');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

lab = xlabel(H.Ax(5,1),'Topography');
pos = get(lab,'position');
pos(1) = 1.05;
set(lab,'position',pos);
set(lab,'units','pixels');

ylabel(H.Ax(1,1),'High')
ylabel(H.Ax(2,1),'Low')
ylabel(H.Ax(3,1),'High')
ylabel(H.Ax(4,1),'Low')
ylabel(H.Ax(5,1),'High')
ylabel(H.Ax(6,1),'Low')

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
xlabel(H.Cbar1,['\Delta till friction angle (' char(176) ')']);
set(H.Cbar1,'ticks',-8:4:8)

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
xlabel(H.Cbar2,'\Delta surface elevation (m)');

%% Plot results
  
% Construct target grounding line contour
H.tempfig = figure;
H.tempax  = axes('parent',H.tempfig);
C_GL_target = contour('parent',H.tempax,'xdata',target.y,'ydata',target.x,'zdata',target.TAF,'levellist',0);
close(H.tempfig);

for i = 1:nay
  
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
      cdata = R.dphi_fric';
      adata = zeros(size(cdata));
      adata( R.TAF'>0) = 1;
      image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled','alphadata',adata);
    end
  
    % Right column: surface elevation
    if j==2
      ax = H.Ax(i,2);
      R  = results(i);
      cdata = R.dHs';
      image('parent',ax,'xdata',R.x,'ydata',R.y,'cdata',cdata,'cdatamapping','scaled');
    end

    % Add "perfect" grounding line contour
    C = C_GL_target;
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
  