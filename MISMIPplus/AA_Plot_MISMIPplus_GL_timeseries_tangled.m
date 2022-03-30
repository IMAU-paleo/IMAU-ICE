clc
clear all
close all

resolutions     = {'2km','5km'};
sliding_laws    = {'Weertman','Tsai2015','Schoof2005'};
stress_balances = {'DIVA','hybrid'};
subgrid_schemes = {'FCMP','PMP','NMP'};

IMAU_ICE.nruns   = 0;
IMAU_ICE.results = [];

for i = 1:length(resolutions    ); IMAU_ICE.(['ens_' resolutions{     i}]).list = []; end
for i = 1:length(sliding_laws   ); IMAU_ICE.(['ens_' sliding_laws{    i}]).list = []; end
for i = 1:length(stress_balances); IMAU_ICE.(['ens_' stress_balances{ i}]).list = []; end
for i = 1:length(subgrid_schemes); IMAU_ICE.(['ens_' subgrid_schemes{ i}]).list = []; end

henk = dir;
for i = 1:length(henk)
  
  if henk(i).isdir && contains(henk(i).name,'MISMIPplus_ice1r_')
    
    IMAU_ICE.nruns = IMAU_ICE.nruns + 1;
    
    % part 1
    filename = [henk(i).name '/MISMIPplus_output.nc'];
    
    time1 = ncread(filename,'time');
    xGL1  = ncread(filename,'xGL');
    ny  = size(xGL1,2);
    jmid = ceil(ny/2);
    xGL1 = xGL1(:,jmid);
    
    % part 2
    filename = ['MISMIPplus_ice1rr' henk(i).name(17:end) '/MISMIPplus_output.nc'];
    
    time2 = ncread(filename,'time');
    xGL2  = ncread(filename,'xGL');
    ny  = size(xGL2,2);
    jmid = ceil(ny/2);
    xGL2 = xGL2(:,jmid);
    
    % combine
    time = [time1(1:end-1);time2];
    xGL  = [xGL1( 1:end-1);xGL2 ];
    
    results.run      = henk(i).name;
    results.time     = time;
    results.xGL      = xGL;
    
    IMAU_ICE.results{end+1} = results;
    
    %% Find out which run this is
    
    % Resolution
    resolution = [];
    for j = 1: length(resolutions)
      if contains(filename,['_' resolutions{j} '_']) || contains(filename,['_' resolutions{j} '/'])
        resolution = resolutions{j};
      end
    end
    if isempty(resolution); error(['Didnt recognise resolution for run "' henk(i).name '"!']); end
    
    % Sliding law
    sliding_law = [];
    for j = 1: length(sliding_laws)
      if contains(filename,['_' sliding_laws{j} '_']) || contains(filename,['_' sliding_laws{j} '/'])
        sliding_law = sliding_laws{j};
      end
    end
    if isempty(sliding_law); error(['Didnt recognise sliding law for run "' henk(i).name '"!']); end
    
    % Stress balance
    stress_balance = [];
    for j = 1: length(stress_balances)
      if contains(filename,['_' stress_balances{j} '_']) || contains(filename,['_' stress_balances{j} '/'])
        stress_balance = stress_balances{j};
      end
    end
    if isempty(stress_balance); error(['Didnt recognise stress balance for run "' henk(i).name '"!']); end
    
    % Subgrid scheme
    subgrid_scheme = [];
    for j = 1: length(subgrid_schemes)
      if contains(filename,['_' subgrid_schemes{j} '_']) || contains(filename,['_' subgrid_schemes{j} '/'])
        subgrid_scheme = subgrid_schemes{j};
      end
    end
    if isempty(subgrid_scheme); error(['Didnt recognise subgrid scheme for run "' henk(i).name '"!']); end
    
    % Add to lists
    IMAU_ICE.(['ens_' resolution    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' sliding_law   ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' stress_balance]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' subgrid_scheme]).list(end+1) = IMAU_ICE.nruns;
    
  end
  
end

%% Create plumes per resolution, sliding law, stress balance, and subgrid scheme
for xi = 1:length(resolutions)
  
  list    = IMAU_ICE.(['ens_' resolutions{xi}]).list;
  n       = length(list);
  
  time    = IMAU_ICE.results{list(1)}.time;
  xGL_av  = zeros(size(time));
  xGL_max = zeros(size(time)) - Inf;
  xGL_min = zeros(size(time)) + Inf;
  
  for i = 1:n
    j = list(i);
    time2 = IMAU_ICE.results{j}.time;
    if any(time2~=time); error('Not all runs have output on the same time frames!'); end
    
    xGL = IMAU_ICE.results{j}.xGL;
    
    xGL_av  = xGL_av + xGL / n;
    xGL_min = min(xGL_min,xGL);
    xGL_max = max(xGL_max,xGL);
  end
  
  IMAU_ICE.(['ens_' resolutions{xi}]).time    = time;
  IMAU_ICE.(['ens_' resolutions{xi}]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' resolutions{xi}]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' resolutions{xi}]).xGL_max = xGL_max;
end
for xi = 1:length(sliding_laws)
  
  list    = IMAU_ICE.(['ens_' sliding_laws{xi}]).list;
  n       = length(list);
  
  time    = IMAU_ICE.results{list(1)}.time;
  xGL_av  = zeros(size(time));
  xGL_max = zeros(size(time)) - Inf;
  xGL_min = zeros(size(time)) + Inf;
  
  for i = 1:n
    j = list(i);
    time2 = IMAU_ICE.results{j}.time;
    if any(time2~=time); error('Not all runs have output on the same time frames!'); end
    
    xGL = IMAU_ICE.results{j}.xGL;
    
    xGL_av  = xGL_av + xGL / n;
    xGL_min = min(xGL_min,xGL);
    xGL_max = max(xGL_max,xGL);
  end
  
  IMAU_ICE.(['ens_' sliding_laws{xi}]).time    = time;
  IMAU_ICE.(['ens_' sliding_laws{xi}]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' sliding_laws{xi}]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' sliding_laws{xi}]).xGL_max = xGL_max;
end
for xi = 1:length(stress_balances)
  
  list    = IMAU_ICE.(['ens_' stress_balances{xi}]).list;
  n       = length(list);
  
  time    = IMAU_ICE.results{list(1)}.time;
  xGL_av  = zeros(size(time));
  xGL_max = zeros(size(time)) - Inf;
  xGL_min = zeros(size(time)) + Inf;
  
  for i = 1:n
    j = list(i);
    time2 = IMAU_ICE.results{j}.time;
    if any(time2~=time); error('Not all runs have output on the same time frames!'); end
    
    xGL = IMAU_ICE.results{j}.xGL;
    
    xGL_av  = xGL_av + xGL / n;
    xGL_min = min(xGL_min,xGL);
    xGL_max = max(xGL_max,xGL);
  end
  
  IMAU_ICE.(['ens_' stress_balances{xi}]).time    = time;
  IMAU_ICE.(['ens_' stress_balances{xi}]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' stress_balances{xi}]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' stress_balances{xi}]).xGL_max = xGL_max;
end
for xi = 1:length(subgrid_schemes)
  
  list    = IMAU_ICE.(['ens_' subgrid_schemes{xi}]).list;
  n       = length(list);
  
  time    = IMAU_ICE.results{list(1)}.time;
  xGL_av  = zeros(size(time));
  xGL_max = zeros(size(time)) - Inf;
  xGL_min = zeros(size(time)) + Inf;
  
  for i = 1:n
    j = list(i);
    time2 = IMAU_ICE.results{j}.time;
    if any(time2~=time); error('Not all runs have output on the same time frames!'); end
    
    xGL = IMAU_ICE.results{j}.xGL;
    
    xGL_av  = xGL_av + xGL / n;
    xGL_min = min(xGL_min,xGL);
    xGL_max = max(xGL_max,xGL);
  end
  
  IMAU_ICE.(['ens_' subgrid_schemes{xi}]).time    = time;
  IMAU_ICE.(['ens_' subgrid_schemes{xi}]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' subgrid_schemes{xi}]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' subgrid_schemes{xi}]).xGL_max = xGL_max;
end

%% Set up GUI

wa          = [    400,    50,   400,    50   ];
margins_hor = [100,    15 ,    15,    15,   15];
ha          = [    300,   300    ];
margins_ver = [80,     40,     40];

nax = length(wa);
nay = length(ha);

wf = sum(margins_hor) + sum(wa);
hf = sum(margins_ver) + sum(ha);

H.Fig = figure('position',[100,100,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    jp = nay+1-j;
    x = sum(margins_hor(1:i )) + sum(wa(1:i -1));
    y = sum(margins_ver(1:jp)) + sum(ha(1:jp-1));
    H.Ax(j,i) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa(i),ha(jp)],'fontsize',24);
  end
end

xlim  = [0,200];
xtick = 0:50:200;
ylim  = [300,475];
ytick = 300:25:475;

set(H.Ax(1,1),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(1,3),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(2,1),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(2,3),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');

set(H.Ax(1,1),'xticklabels',{});
set(H.Ax(1,3),'xticklabels',{});
set(H.Ax(1,3),'yticklabels',{});
set(H.Ax(2,3),'yticklabels',{});

title(H.Ax(1,1),'A');
title(H.Ax(1,3),'B');
title(H.Ax(2,1),'C');
title(H.Ax(2,3),'D');

H.Ax(1,2).XAxis.Visible = 'off';
H.Ax(1,4).XAxis.Visible = 'off';
H.Ax(2,2).XAxis.Visible = 'off';
H.Ax(2,4).XAxis.Visible = 'off';
H.Ax(1,2).YAxis.Visible = 'off';
H.Ax(1,4).YAxis.Visible = 'off';
H.Ax(2,2).YAxis.Visible = 'off';
H.Ax(2,4).YAxis.Visible = 'off';
set(H.Ax(1,2),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,4),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(2,2),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(2,4),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')

xlabel(H.Ax(2,1),'Time (yr)')
xlabel(H.Ax(2,3),'Time (yr)')
ylabel(H.Ax(1,1),'x_{GL} (km)')
ylabel(H.Ax(2,1),'x_{GL} (km)')

%% Plot plumes

% Resolution
plumes = [];
plumes{1} = IMAU_ICE.ens_5km;
plumes{2} = IMAU_ICE.ens_2km;
plot_plumes(H.Ax(1,1),H.Ax(1,2),plumes)
legend(H.Ax(1,1),'5 km','2 km','location','northeast')

% Sliding law
plumes = [];
plumes{1} = IMAU_ICE.ens_Weertman;
plumes{2} = IMAU_ICE.ens_Tsai2015;
plumes{3} = IMAU_ICE.ens_Schoof2005;
plot_plumes(H.Ax(1,3),H.Ax(1,4),plumes)
legend(H.Ax(1,3),'Weertman','Tsai2015','Schoof2005','location','northeast')

% Stress balance
plumes = [];
plumes{1} = IMAU_ICE.ens_DIVA;
plumes{2} = IMAU_ICE.ens_hybrid;
plot_plumes(H.Ax(2,1),H.Ax(2,2),plumes)
legend(H.Ax(2,1),'DIVA','hybrid','location','northeast')

% Subgrid scheme
plumes = [];
plumes{1} = IMAU_ICE.ens_NMP;
plumes{2} = IMAU_ICE.ens_FCMP;
plumes{3} = IMAU_ICE.ens_PMP;
plot_plumes(H.Ax(2,3),H.Ax(2,4),plumes)
legend(H.Ax(2,3),'NMP','FCMP','PMP','location','northeast')

function plot_plumes(ax1,ax2,plumes)

  n = length(plumes);
  colors = lines(n);
  
  % Empty line objects for legend
  for i = 1:length(plumes)
    line('parent',ax1,'xdata',[],'ydata',[],'color',colors(i,:),'linewidth',3);
  end

  for i = 1:length(plumes)
    
    % Plot plume in big axes
    time    = plumes{i}.time;
    xGL_av  = plumes{i}.xGL_av;
    xGL_min = plumes{i}.xGL_min;
    xGL_max = plumes{i}.xGL_max;
    
    xdata_lo =        time;
    xdata_hi = flipud(time);
    ydata_lo =        xGL_min;
    ydata_hi = flipud(xGL_max);
    xdata = [xdata_lo;xdata_hi];
    ydata = [ydata_lo;ydata_hi];
    
    patch('parent',ax1,'xdata',xdata,'ydata',ydata/1e3,'facecolor',colors(i,:),'facealpha',0.25,'edgecolor','none');
    line('parent',ax1,'xdata',time,'ydata',xGL_av/1e3,'color',colors(i,:),'linewidth',3);
    
    % Plot plume end in small axes
    w = 1/n;
    m = -w/2 + i*w;
    
    xl = m - w*0.4;
    xr = m + w*0.4;
    
    yl = xGL_min(end)/1e3;
    yu = xGL_max(end)/1e3;
    
    patch('parent',ax2,'xdata',[xl,xr,xr,xl],'ydata',[yl,yl,yu,yu],'facecolor',colors(i,:),'facealpha',0.25,'edgecolor','none');
    line('parent',ax2,'xdata',[xl,xr],'ydata',[0,0]+xGL_av(end)/1e3,'color',colors(i,:),'linewidth',3);
  end
  
end