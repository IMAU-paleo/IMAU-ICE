clc
clear all
close all

clim_Hi  = [0,3500];
clim_phi = [0,10];

%% Read data
foldernames = {...
  'BIVMIP_A_perfect_10km',...
  'BIVMIP_A_inv_10km_perfect',...
  'BIVMIP_A_inv_10km_visc_lo',...
  'BIVMIP_A_inv_10km_visc_hi',...
  'BIVMIP_A_inv_10km_SMB_lo',...
  'BIVMIP_A_inv_10km_SMB_hi',...
  'BIVMIP_A_inv_10km_p_lo',...
  'BIVMIP_A_inv_10km_p_hi',...
  'BIVMIP_A_inv_10km_ut_hi',...
  'BIVMIP_A_inv_10km_ut_hi'};

for fi = 1: length(foldernames)
  
  filename_restart       = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields   = [foldernames{fi} '/help_fields_ANT.nc'];
  
  results(fi).x         = ncread(filename_restart    ,'x');
  results(fi).y         = ncread(filename_restart    ,'y');
  results(fi).time      = ncread(filename_restart    ,'time');
  ti = length(results(fi).time);
  
  results(fi).Hi        = ncread(filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).phi_fric  = ncread(filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
end

%% Set up GUI
wa = 180;
ha = 180;

margins_hor = [75,25,25,25,25,25,120];
margins_ver = [25,25,25,25,50];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

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
      'fontsize',24,'xlim',[-700e3,700e3],'xtick',[],'ylim',[-700e3,700e3],'ytick',[]);
    
    if (i==1)
      set(H.Ax(i,j),'xaxislocation','top')
    end
    
    if (i==1 || i==2)
      colormap(H.Ax(i,j),parula(16));
      set(H.Ax(i,j),'clim',clim_phi);
    else
      colormap(H.Ax(i,j),itmap(16));
      set(H.Ax(i,j),'clim',clim_Hi);
    end
  end
end

ylabel(H.Ax(1,1),['Till friction angle (' char(176) ')']);
ylabel(H.Ax(3,1),'Ice thickness (m)');

xlabel(H.Ax(1,1),'Analytical');
xlabel(H.Ax(1,2),'Perfect');
xlabel(H.Ax(1,3),'Viscosity');
xlabel(H.Ax(1,4),'SMB');
xlabel(H.Ax(1,5),'Z-I p');
xlabel(H.Ax(1,6),'Z-I u\_t');

%% Colorbars

% phi
pos1 = get(H.Ax(1,nax),'position');
pos2 = get(H.Ax(2,nax),'position');

x    = pos1(1) + pos1(3) + 20;
ymin = pos2(2);
ymax = pos1(2) + pos1(3);

H.Axcbar1 = axes('parent',H.Fig,'units','pixels','position',[x,ymin,100,ymax-ymin],...
  'xtick',[],'ytick',[],'box','off','fontsize',24,'clim',clim_phi);
H.Axcbar1.XAxis.Visible = 'off';
H.Axcbar1.YAxis.Visible = 'off';
colormap(H.Axcbar1,parula(16))
H.Cbar1 = colorbar(H.Axcbar1,'location','west');
ylabel(H.Cbar1,['Till friction angle (' char(176) ')'])

% Hi
pos1 = get(H.Ax(3,nax),'position');
pos2 = get(H.Ax(4,nax),'position');

x    = pos1(1) + pos1(3) + 20;
ymin = pos2(2);
ymax = pos1(2) + pos1(3);

H.Axcbar2 = axes('parent',H.Fig,'units','pixels','position',[x,ymin,100,ymax-ymin],...
  'xtick',[],'ytick',[],'box','off','fontsize',24,'clim',clim_Hi);
H.Axcbar2.XAxis.Visible = 'off';
H.Axcbar2.YAxis.Visible = 'off';
colormap(H.Axcbar2,itmap(16))
H.Cbar2 = colorbar(H.Axcbar2,'location','west');
ylabel(H.Cbar2,'Ice thickness (m)')

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
  
% % Mark ice-free area as uncertain
% for j = 1:nax
%   image('parent',H.Ax(1,j),'xdata',results(j).x,'ydata',results(j).y,'cdata',ones(length(results(j).x),length(results(j).y),3),...
%     'alphadata',double(results(j).Hi'<0.1)*0.5,'alphadatamapping','none');
% end
