clc
clear all
close all

do_plot_errors = false;

%% Read data
foldernames = {...
  'BIVMIP_B_perfect_10km',...
  'BIVMIP_B_inv_10km_perfect',...
  'BIVMIP_B_inv_10km_visc_lo',...
  'BIVMIP_B_inv_10km_visc_hi',...
  'BIVMIP_B_inv_10km_SMB_lo',...
  'BIVMIP_B_inv_10km_SMB_hi',...
  'BIVMIP_B_inv_10km_ut_hi',...
  'BIVMIP_B_inv_10km_ut_lo',...
  'BIVMIP_B_inv_10km_p_lo',...
  'BIVMIP_B_inv_10km_p_hi'};

clim_Hi  = [0,3000];
clim_phi = [0,7];

phi_map  = parula(256);
Hi_map   = itmap(16);

clim_dHi = [-100,100];
clim_dphi = [-5,5];

dphi_map = flipud(lbmap(64,'redblue'));
dHi_map  = flipud(lbmap(64,'redblue'));

for fi = 1: length(foldernames)
  
  filename_restart       = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields   = [foldernames{fi} '/help_fields_ANT.nc'];
  
  results(fi).x         = ncread(filename_restart    ,'x');
  results(fi).y         = ncread(filename_restart    ,'y');
  results(fi).time      = ncread(filename_restart    ,'time');
  ti = length(results(fi).time);
  
  results(fi).Hi        = ncread(filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).phi_fric  = ncread(filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  
  if fi > 1 && do_plot_errors
    results( fi).Hi       = results( fi).Hi       - results( 1).Hi;
    results( fi).phi_fric = results( fi).phi_fric - results( 1).phi_fric;
  end
end

%% Set up GUI
wa = 180;
ha = 180;

margins_hor = [75,5,5,5,5,5,120];
margins_ver = [25,5,5,5,50];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

xlim = [-575,575] * 1e3;
ylim = [-600,550] * 1e3;

H.Fig = figure('position',[100,200,wf,hf],'color','w');
H.Ax  = zeros( nay,nax);
H.Axa = zeros( nay,nax);
for i = 1: nay
  for j = 1: nax
    
    if (i == 2 && j == 1); continue; end
    if (i == 2 && j == 2); continue; end
    
    if (i == 4 && j == 1); continue; end
    if (i == 4 && j == 2); continue; end
    
    x = sum(margins_hor(1:j )) + (j -1)*wa;
    ip = nay+1-i;
    y = sum(margins_ver(1:ip)) + (ip-1)*ha;
    H.Ax( i,j) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],...
      'fontsize',24,'xlim',xlim,'xtick',[],'ylim',ylim,'ytick',[]);
    
    if (i==1)
      set(H.Ax(i,j),'xaxislocation','top')
    end
    
    if (i==1 || i==2)
      if j == 1 || ~do_plot_errors
        colormap(H.Ax(i,j),phi_map);
        set(H.Ax(i,j),'clim',clim_phi);
      else
        colormap(H.Ax(i,j),dphi_map);
        set(H.Ax(i,j),'clim',clim_dphi);
      end
    else
      if j == 1 || ~do_plot_errors
        colormap(H.Ax(i,j),Hi_map);
        set(H.Ax(i,j),'clim',clim_Hi);
      else
        colormap(H.Ax(i,j),dHi_map);
        set(H.Ax(i,j),'clim',clim_dHi);
      end
    end
  end
end

ylabel(H.Ax(1,1),['Till friction angle (' char(176) ')']);
ylabel(H.Ax(3,1),'Ice thickness (m)');

xlabel(H.Ax(1,1),'Analytical');
xlabel(H.Ax(1,2),'Perfect');
xlabel(H.Ax(1,3),'Viscosity');
xlabel(H.Ax(1,4),'SMB');
xlabel(H.Ax(1,5),'Z-I u\_t');
xlabel(H.Ax(1,6),'Z-I p');

%% Colorbars

% phi
pos1 = get(H.Ax(1,nax),'position');
pos2 = get(H.Ax(2,nax),'position');

x    = pos1(1) + pos1(3) + 20;
ymin = pos2(2);
ymax = pos1(2) + pos1(3);

H.Axcbar1 = axes('parent',H.Fig,'units','pixels','position',[x,ymin,100,ymax-ymin],...
  'xtick',[],'ytick',[],'box','off','fontsize',24);
H.Axcbar1.XAxis.Visible = 'off';
H.Axcbar1.YAxis.Visible = 'off';
if ~do_plot_errors
  set(H.Axcbar1,'clim',clim_phi)
  colormap(H.Axcbar1,phi_map)
else
  set(H.Axcbar1,'clim',clim_dphi)
  colormap(H.Axcbar1,dphi_map)
end
H.Cbar1 = colorbar(H.Axcbar1,'location','west');
if ~do_plot_errors
  ylabel(H.Cbar1,['Till friction angle (' char(176) ')'])
else
  ylabel(H.Cbar1,['\Delta Till friction angle (' char(176) ')'])
end

% Hi
pos1 = get(H.Ax(3,nax),'position');
pos2 = get(H.Ax(4,nax),'position');

x    = pos1(1) + pos1(3) + 20;
ymin = pos2(2);
ymax = pos1(2) + pos1(3);

H.Axcbar2 = axes('parent',H.Fig,'units','pixels','position',[x,ymin,100,ymax-ymin],...
  'xtick',[],'ytick',[],'box','off','fontsize',24);
H.Axcbar2.XAxis.Visible = 'off';
H.Axcbar2.YAxis.Visible = 'off';
if ~do_plot_errors
  set(H.Axcbar2,'clim',clim_Hi)
  colormap(H.Axcbar2,Hi_map)
else
  set(H.Axcbar2,'clim',clim_dHi)
  colormap(H.Axcbar2,dHi_map)
end
H.Cbar2 = colorbar(H.Axcbar2,'location','west');
if ~do_plot_errors
  ylabel(H.Cbar2,'Ice thickness (m)')
else
  ylabel(H.Cbar2,'\Delta Ice thickness (m)')
end

%% Plot data

% "Perfect" forward solution
image('parent',H.Ax(1,1),'xdata',results(1).x,'ydata',results(1).y,'cdata',results(1).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,1),'xdata',results(1).x,'ydata',results(1).y,'cdata',results(1).Hi','cdatamapping','scaled');

% "Perfect" inversion
image('parent',H.Ax(1,2),'xdata',results(2).x,'ydata',results(2).y,'cdata',results(2).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,2),'xdata',results(2).x,'ydata',results(2).y,'cdata',results(2).Hi','cdatamapping','scaled');

% "Perturbed" inversions
j = 3;
n = 3;
image('parent',H.Ax(1,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');
n = 4;
image('parent',H.Ax(2,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(4,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');

j = 4;
n = 5;
image('parent',H.Ax(1,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');
n = 6;
image('parent',H.Ax(2,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(4,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');

j = 5;
n = 7;
image('parent',H.Ax(1,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');
n = 8;
image('parent',H.Ax(2,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(4,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');

j = 6;
n = 9;
image('parent',H.Ax(1,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(3,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');
n = 10;
image('parent',H.Ax(2,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).phi_fric','cdatamapping','scaled');
image('parent',H.Ax(4,j),'xdata',results(n).x,'ydata',results(n).y,'cdata',results(n).Hi','cdatamapping','scaled');
  
% Mark ice-free area as uncertain
alpha_noicemask = 1.0;
for j = 1:nax
  if (H.Ax(1,j)>0)
    image('parent',H.Ax(1,j),'xdata',results(1).x,'ydata',results(1).y,'cdata',ones(length(results(1).x),length(results(1).y),3),...
      'alphadata',double(results(1).Hi'<0.1)*alpha_noicemask,'alphadatamapping','none');
  end
  if (H.Ax(2,j)>0)
    image('parent',H.Ax(2,j),'xdata',results(1).x,'ydata',results(1).y,'cdata',ones(length(results(1).x),length(results(1).y),3),...
      'alphadata',double(results(1).Hi'<0.1)*alpha_noicemask,'alphadatamapping','none');
  end
end
