clc
clear all
close all

experiments     = {'IceOcean0','IceOcean1rr','IceOcean1ra','IceOcean2rr','IceOcean2ra'};
resolutions     = {'2km','5km'};
subgrid_schemes = {'FCMP','PMP','NMP'};
melt_models     = {'lin','quad','Mplus','plume','PICO','PICOP'};

IMAU_ICE.nruns   = 0;
IMAU_ICE.results = [];

for i = 1:length(experiments    ); IMAU_ICE.(['ens_' experiments{     i}]).list = []; end
for i = 1:length(resolutions    ); IMAU_ICE.(['ens_' resolutions{     i}]).list = []; end
for i = 1:length(subgrid_schemes); IMAU_ICE.(['ens_' subgrid_schemes{ i}]).list = []; end
for i = 1:length(melt_models    ); IMAU_ICE.(['ens_' melt_models{     i}]).list = []; end

% Plumes per [experiment + melt model]
for i = 1:length(experiments)
  for j = 1:length(melt_models)
    IMAU_ICE.(['ens_' experiments{i} '_' melt_models{j}]).list = [];
  end
end

henk = dir;
for i = 1:length(henk)
  
  if henk(i).isdir && contains(henk(i).name,'MISOMIP1_IceOcean')
    
    filename = [henk(i).name '/MISMIPplus_output.nc'];
    
    % DENK DROM
    if ~exist(filename,'file'); continue; end
    
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
    
    % Experiment
    experiment = [];
    for j = 1: length(experiments)
      if contains(filename,['_' experiments{j} '_']) || contains(filename,['_' experiments{j} '/'])
        experiment = experiments{j};
      end
    end
    if isempty(experiment); error(['Didnt recognise experiment for run "' henk(i).name '"!']); end
    
    % Resolution
    resolution = [];
    for j = 1: length(resolutions)
      if contains(filename,['_' resolutions{j} '_']) || contains(filename,['_' resolutions{j} '/'])
        resolution = resolutions{j};
      end
    end
    if isempty(resolution); error(['Didnt recognise resolution for run "' henk(i).name '"!']); end
    
    % Subgrid scheme
    subgrid_scheme = [];
    for j = 1: length(subgrid_schemes)
      if contains(filename,['_' subgrid_schemes{j} '_']) || contains(filename,['_' subgrid_schemes{j} '/'])
        subgrid_scheme = subgrid_schemes{j};
      end
    end
    if isempty(subgrid_scheme); error(['Didnt recognise subgrid scheme for run "' henk(i).name '"!']); end
    
    % Melt model
    melt_model = [];
    for j = 1: length(melt_models)
      if contains(filename,['_' melt_models{j} '_']) || contains(filename,['_' melt_models{j} '/'])
        melt_model = melt_models{j};
      end
    end
    if isempty(melt_model); error(['Didnt recognise subgrid scheme for run "' henk(i).name '"!']); end
    
    % Add to lists
    IMAU_ICE.(['ens_' experiment    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' resolution    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' subgrid_scheme]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' melt_model    ]).list(end+1) = IMAU_ICE.nruns;
    
    % Plumes per [experiment + melt model]
    IMAU_ICE.(['ens_' experiment '_' melt_model]).list(end+1) = IMAU_ICE.nruns;
    
  end
  
end

% Calculate plumes per [experiment + melt model]
for xi = 1:length(experiments)
  experiment = experiments{xi};
  for mi = 1:length(melt_models)
    melt_model = melt_models{mi};
  
    list    = IMAU_ICE.(['ens_' experiment '_' melt_model]).list;
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

    IMAU_ICE.(['ens_' experiment '_' melt_model]).time    = time;
    IMAU_ICE.(['ens_' experiment '_' melt_model]).xGL_av  = xGL_av;
    IMAU_ICE.(['ens_' experiment '_' melt_model]).xGL_min = xGL_min;
    IMAU_ICE.(['ens_' experiment '_' melt_model]).xGL_max = xGL_max;
  
  end
end

%% Set up GUI

wa          = [    400,    40,   40,  40,   400,    40,  40,  40   ];
margins_hor = [100,    15 ,    0,   0,   15,    15,    0,   0,   15];
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
set(H.Ax(1,5),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');

set(H.Ax(1,5),'yticklabels',{});

title(H.Ax(1,1),'A');
title(H.Ax(1,5),'B');

H.Ax(1,2).XAxis.Visible = 'off';
H.Ax(1,3).XAxis.Visible = 'off';
H.Ax(1,4).XAxis.Visible = 'off';
H.Ax(1,6).XAxis.Visible = 'off';
H.Ax(1,7).XAxis.Visible = 'off';
H.Ax(1,8).XAxis.Visible = 'off';

H.Ax(1,2).YAxis.Visible = 'off';
H.Ax(1,3).YAxis.Visible = 'off';
H.Ax(1,4).YAxis.Visible = 'off';
H.Ax(1,6).YAxis.Visible = 'off';
H.Ax(1,7).YAxis.Visible = 'off';
H.Ax(1,8).YAxis.Visible = 'off';

set(H.Ax(1,2),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,3),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,4),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,6),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,7),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')
set(H.Ax(1,8),'ylim',ylim,'ytick',ytick,'ygrid','on','xlim',[0,1],'xtick',[0,1],'xgrid','on')

xlabel(H.Ax(1,1),'Time (yr)')
xlabel(H.Ax(1,5),'Time (yr)')
ylabel(H.Ax(1,1),'x_{GL} (km)')

text(H.Ax(1,2),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'0' ,'fontsize',24)
text(H.Ax(1,3),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'ra','fontsize',24)
text(H.Ax(1,4),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'rr','fontsize',24)
text(H.Ax(1,6),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'0' ,'fontsize',24)
text(H.Ax(1,7),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'ra','fontsize',24)
text(H.Ax(1,8),0.3,ylim(2)+(ylim(2)-ylim(1))*0.05,'rr','fontsize',24)

%% Plot

% IceOcean0
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean0_lin;
plumes{2} = IMAU_ICE.ens_IceOcean0_quad;
plumes{3} = IMAU_ICE.ens_IceOcean0_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean0_plume;
plumes{5} = IMAU_ICE.ens_IceOcean0_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean0_PICOP;
plot_plumes(H.Ax(1,1),H.Ax(1,2),plumes)

% IceOcean1ra
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean1ra_lin;
plumes{2} = IMAU_ICE.ens_IceOcean1ra_quad;
plumes{3} = IMAU_ICE.ens_IceOcean1ra_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean1ra_plume;
plumes{5} = IMAU_ICE.ens_IceOcean1ra_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean1ra_PICOP;
for i = 1:length(plumes)
  m = plumes{i}.time<100;
  plumes{i}.time(    m) = [];
  plumes{i}.xGL_av(  m) = [];
  plumes{i}.xGL_min( m) = [];
  plumes{i}.xGL_max( m) = [];
end
plot_plumes(H.Ax(1,1),H.Ax(1,3),plumes)

% IceOcean1rr
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean1rr_lin;
plumes{2} = IMAU_ICE.ens_IceOcean1rr_quad;
plumes{3} = IMAU_ICE.ens_IceOcean1rr_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean1rr_plume;
plumes{5} = IMAU_ICE.ens_IceOcean1rr_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean1rr_PICOP;
plot_plumes(H.Ax(1,1),H.Ax(1,4),plumes)

legend(H.Ax(1,1),'Linear','Quadratic','M+','Plume','PICO','PICOP','location','southwest')

% IceOcean0
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean0_lin;
plumes{2} = IMAU_ICE.ens_IceOcean0_quad;
plumes{3} = IMAU_ICE.ens_IceOcean0_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean0_plume;
plumes{5} = IMAU_ICE.ens_IceOcean0_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean0_PICOP;
plot_plumes(H.Ax(1,5),H.Ax(1,6),plumes)

% IceOcean2ra
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean2ra_lin;
plumes{2} = IMAU_ICE.ens_IceOcean2ra_quad;
plumes{3} = IMAU_ICE.ens_IceOcean2ra_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean2ra_plume;
plumes{5} = IMAU_ICE.ens_IceOcean2ra_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean2ra_PICOP;
for i = 1:length(plumes)
  m = plumes{i}.time<100;
  plumes{i}.time(    m) = [];
  plumes{i}.xGL_av(  m) = [];
  plumes{i}.xGL_min( m) = [];
  plumes{i}.xGL_max( m) = [];
end
plot_plumes(H.Ax(1,5),H.Ax(1,7),plumes)

% IceOcean2rr
plumes = [];
plumes{1} = IMAU_ICE.ens_IceOcean2rr_lin;
plumes{2} = IMAU_ICE.ens_IceOcean2rr_quad;
plumes{3} = IMAU_ICE.ens_IceOcean2rr_Mplus;
plumes{4} = IMAU_ICE.ens_IceOcean2rr_plume;
plumes{5} = IMAU_ICE.ens_IceOcean2rr_PICO;
plumes{6} = IMAU_ICE.ens_IceOcean2rr_PICOP;
plot_plumes(H.Ax(1,5),H.Ax(1,8),plumes)

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