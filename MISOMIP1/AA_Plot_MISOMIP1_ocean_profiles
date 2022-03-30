clc
clear all
close all

% Cold
cold.T0   = -1.9;
cold.Tbot = -1.9;
cold.S0   = 33.8;
cold.Sbot = 34.55;

% Warm
warm.T0   = -1.9;
warm.Tbot =  1.0;
warm.S0   = 33.8;
warm.Sbot = 34.7;

depth_max = 720;
depth = linspace(0,1000,1000);

cold.T = cold.T0 + max(0,min(1, depth/depth_max)) * (cold.Tbot - cold.T0);
warm.T = warm.T0 + max(0,min(1, depth/depth_max)) * (warm.Tbot - warm.T0);
cold.S = cold.S0 + max(0,min(1, depth/depth_max)) * (cold.Sbot - cold.S0);
warm.S = warm.S0 + max(0,min(1, depth/depth_max)) * (warm.Sbot - warm.S0);

%% Plot
wa = 300;
ha = 300;

margins_hor = [100,40,25];
margins_ver = [75,25];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = zeros( nay,nax);
H.Axa = zeros( nay,nax);
for i = 1: nay
  for j = 1: nax
    x = sum(margins_hor(1:j )) + (j -1)*wa;
    ip = nay+1-i;
    y = sum(margins_ver(1:ip)) + (ip-1)*ha;
    H.Ax( i,j) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24,...
      'xgrid','on','ygrid','on');
  end
end

set(H.Ax(1,1),'xlim',[ -2.1, 1.1],'ylim',[0,1000],'ydir','reverse');
set(H.Ax(1,2),'xlim',[ 33.5,35.0],'ylim',[0,1000],'ydir','reverse','yticklabels',[]);

ylabel(H.Ax(1,1),'Depth (m)');
xlabel(H.Ax(1,1),['Temperature (' char(176) 'C)']);
xlabel(H.Ax(1,2),'Salinity (PSU)');

% Empty line objects for legend
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'linewidth',3,'color','b');
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'linewidth',3,'color','r');

% Plot data
xlim = get(H.Ax(1,1),'xlim');
line('parent',H.Ax(1,1),'xdata',xlim,'ydata',[0,0]+depth_max,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',cold.T,'ydata',depth,'linewidth',3,'color','b');
line('parent',H.Ax(1,1),'xdata',warm.T,'ydata',depth,'linewidth',3,'color','r');

xlim = get(H.Ax(1,2),'xlim');
line('parent',H.Ax(1,2),'xdata',xlim,'ydata',[0,0]+depth_max,'linestyle','--');
line('parent',H.Ax(1,2),'xdata',cold.S,'ydata',depth,'linewidth',3,'color','b');
line('parent',H.Ax(1,2),'xdata',warm.S,'ydata',depth,'linewidth',3,'color','r');

% Legend
legend(H.Ax(1,1),'COLD','WARM');