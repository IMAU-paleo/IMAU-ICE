clc
clear all
close all

resolutions     = {'2km','5km'};
sliding_laws    = {'Weertman','Tsai2015','Schoof2005'};
melt_models     = {'lin','quad','Mplus','plume','PICO','PICOP'};

IMAU_ICE.nruns   = 0;
IMAU_ICE.results = [];

for i = 1:length(resolutions    ); IMAU_ICE.(['ens_' resolutions{     i}]).list = []; end
for i = 1:length(sliding_laws   ); IMAU_ICE.(['ens_' sliding_laws{    i}]).list = []; end
for i = 1:length(melt_models    ); IMAU_ICE.(['ens_' melt_models{     i}]).list = []; end

henk = dir;
for i = 1:length(henk)
  
  if henk(i).isdir && contains(henk(i).name,'MISOMIP1_IceOcean1rr')
    
    filename = [henk(i).name '/MISMIPplus_output.nc'];
    
    IMAU_ICE.nruns = IMAU_ICE.nruns + 1;
    
    time = ncread(filename,'time');
    xGL  = ncread(filename,'xGL');
    ny   = size(xGL,2);
    jmid = ceil(ny/2);
    xGL  = xGL(:,jmid);
    
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
    
    % Subgrid scheme
    sliding_law = [];
    for j = 1: length(sliding_laws)
      if contains(filename,['_' sliding_laws{j} '_']) || contains(filename,['_' sliding_laws{j} '/'])
        sliding_law = sliding_laws{j};
      end
    end
    if isempty(sliding_law); error(['Didnt recognise subgrid scheme for run "' henk(i).name '"!']); end
    
    % Melt model
    melt_model = [];
    for j = 1: length(melt_models)
      if contains(filename,['_' melt_models{j} '_']) || contains(filename,['_' melt_models{j} '/'])
        melt_model = melt_models{j};
      end
    end
    if isempty(melt_model); error(['Didnt recognise subgrid scheme for run "' henk(i).name '"!']); end
    
    % Add to lists
    IMAU_ICE.(['ens_' resolution    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' sliding_law   ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' melt_model    ]).list(end+1) = IMAU_ICE.nruns;
    
  end
  
end

% Calculate plumes
for xi = 1:length(resolutions)
  resolution = resolutions{xi};
  
  list    = IMAU_ICE.(['ens_' resolution]).list;
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

  IMAU_ICE.(['ens_' resolution]).time    = time;
  IMAU_ICE.(['ens_' resolution]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' resolution]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' resolution]).xGL_max = xGL_max;
  
end
for xi = 1:length(sliding_laws)
  sliding_law = sliding_laws{xi};
  
  list    = IMAU_ICE.(['ens_' sliding_law]).list;
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

  IMAU_ICE.(['ens_' sliding_law]).time    = time;
  IMAU_ICE.(['ens_' sliding_law]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' sliding_law]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' sliding_law]).xGL_max = xGL_max;
  
end
for xi = 1:length(melt_models)
  melt_model = melt_models{xi};
  
  list    = IMAU_ICE.(['ens_' melt_model]).list;
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

  IMAU_ICE.(['ens_' melt_model]).time    = time;
  IMAU_ICE.(['ens_' melt_model]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' melt_model]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' melt_model]).xGL_max = xGL_max;
  
end

%% Set up GUI

wa          = [    400,    40,   400,    40,    400,   40   ];
margins_hor = [100,    15 ,    15,    15,    15,    15,   15];
ha          = [    300   ];
margins_ver = [80,     40];

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
set(H.Ax(1,5),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');

set(H.Ax(1,3),'yticklabels',{});
set(H.Ax(1,5),'yticklabels',{});

title(H.Ax(1,1),'A');
title(H.Ax(1,3),'B');
title(H.Ax(1,5),'C');

H.Ax(1,2).XAxis.Visible = 'off';
H.Ax(1,4).XAxis.Visible = 'off';
H.Ax(1,6).XAxis.Visible = 'off';

H.Ax(1,2).YAxis.Visible = 'off';
H.Ax(1,4).YAxis.Visible = 'off';
H.Ax(1,6).YAxis.Visible = 'off';

set(H.Ax(1,2),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,4),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,6),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')

xlabel(H.Ax(1,1),'Time (yr)')
xlabel(H.Ax(1,3),'Time (yr)')
xlabel(H.Ax(1,5),'Time (yr)')
ylabel(H.Ax(1,1),'x_{GL} (km)')

%% Plot

% Resolutions
plumes = [];
plumes{1} = IMAU_ICE.ens_5km;
plumes{2} = IMAU_ICE.ens_2km;
plot_plumes(H.Ax(1,1),H.Ax(1,2),plumes)
legend(H.Ax(1,1),'5 km','2 km','location','southwest')

% Sliding laws
plumes = [];
plumes{1} = IMAU_ICE.ens_Weertman;
plumes{2} = IMAU_ICE.ens_Tsai2015;
plumes{3} = IMAU_ICE.ens_Schoof2005;
plot_plumes(H.Ax(1,3),H.Ax(1,4),plumes)
legend(H.Ax(1,3),'Weertman','Tsai2015','Schoof2005','location','southwest')

% Melt models
plumes = [];
plumes{1} = IMAU_ICE.ens_lin;
plumes{2} = IMAU_ICE.ens_quad;
plumes{3} = IMAU_ICE.ens_Mplus;
plumes{4} = IMAU_ICE.ens_plume;
plumes{5} = IMAU_ICE.ens_PICO;
plumes{6} = IMAU_ICE.ens_PICOP;
plot_plumes(H.Ax(1,5),H.Ax(1,6),plumes)
legend(H.Ax(1,5),'Linear','Quadratic','M+','Plume','PICO','PICOP','location','southwest')


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
    
    patch('parent',ax2,'xdata',[xl,xr,xr,xl],'ydata',[yl,yl,yu,yu],'facecolor',colors(i,:),'facealpha',0.35,'edgecolor','none');
    line('parent',ax2,'xdata',[xl,xr],'ydata',[0,0]+xGL_av(end)/1e3,'color',colors(i,:),'linewidth',3);
  end
  
end