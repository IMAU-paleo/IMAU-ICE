clc
clear all
close all

%% Collect data from IMAU-ICE

experiments     = {'IceOcean0','IceOcean1ra','IceOcean1rr','IceOcean2ra','IceOcean2rr'};
resolutions     = {'5km','2km'};
sliding_laws    = {'slid1','slid2','slid3'};
melt_models     = {'lin','quad','Mplus','plume','PICO','PICOP'};
subgrid_schemes = {'NMP','FCMP','PMP'};

% MISOMIP1
for xi = 1:length(experiments)
  xp = experiments{xi};
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};
      for mi = 1:length(melt_models)
        melt = melt_models{mi};
      
        foldername = ['MISOMIP1_' xp '_' res '_' slid '_' melt];
        filename = [foldername '/MISMIPplus_output.nc'];

        if exist(filename,'file')
          disp(['Reading IMAU-ICE data for run ' foldername '...']);
          time = ncread(filename,'time');
          xGL  = ncread(filename,'xGL');
          ny  = size(xGL,2);
          jmid = ceil(ny/2);
          xGL = xGL(:,jmid);
        else
          time = [];
          xGL  = [];
        end

        MISOMIP1.(xp).(rres).(slid).(melt).time = time;
        MISOMIP1.(xp).(rres).(slid).(melt).xGL  = xGL;
      
      end
    end
  end
end

% Subgridmelt
for ri = 1:length(resolutions)
  res  = resolutions{ri};
  rres = ['r' resolutions{ri}];
  for mi = 1:length(melt_models)
    melt = melt_models{mi};
    for sbgi = 1:length(subgrid_schemes)
      sbg = subgrid_schemes{sbgi};

      foldername = ['subgridmelt_' res '_' melt '_' sbg];
      filename = [foldername '/MISMIPplus_output.nc'];

      if exist(filename,'file')
        disp(['Reading IMAU-ICE data for run ' foldername '...']);
        time = ncread(filename,'time');
        xGL  = ncread(filename,'xGL');
        ny  = size(xGL,2);
        jmid = ceil(ny/2);
        xGL = xGL(:,jmid);
      else
        time = [];
        xGL  = [];
      end

      subgridmelt.(rres).(melt).(sbg).time = time;
      subgridmelt.(rres).(melt).(sbg).xGL  = xGL;

    end
  end
end

%% Create plumes per experiment per melt model

model.time    = (0:5:200)';
model.xGL_av  = zeros(size(time));
model.xGL_min = zeros(size(time)) + Inf;
model.xGL_max = zeros(size(time)) - Inf;
model.nmod    = 0;

for mi = 1:length(melt_models)
  melt = melt_models{mi};
  
  for xi = 1:length(experiments)
    xp = experiments{xi};
  
    MISOMIP1.(melt).(xp) = model;

    for ri = 1:length(resolutions)
      res  = resolutions{ri};
      rres = ['r' resolutions{ri}];
      for sli = 1:length(sliding_laws)
        slid = sliding_laws{sli};
          
        if ~isempty(MISOMIP1.(xp).(rres).(slid).(melt).time)
          MISOMIP1.(melt).(xp).nmod    =      MISOMIP1.(melt).(xp).nmod + 1;
          MISOMIP1.(melt).(xp).xGL_av  =      MISOMIP1.(melt).(xp).xGL_av + MISOMIP1.(xp).(rres).(slid).(melt).xGL;
          MISOMIP1.(melt).(xp).xGL_min = min( MISOMIP1.(melt).(xp).xGL_min, MISOMIP1.(xp).(rres).(slid).(melt).xGL);
          MISOMIP1.(melt).(xp).xGL_max = max( MISOMIP1.(melt).(xp).xGL_max, MISOMIP1.(xp).(rres).(slid).(melt).xGL);
        end
          
      end
    end
  
    if MISOMIP1.(melt).(xp).nmod > 0
      MISOMIP1.(melt).(xp).xGL_av = MISOMIP1.(melt).(xp).xGL_av / MISOMIP1.(melt).(xp).nmod;
    end
    
  end
end

%% Create plumes for sliding law and sub-grid scheme per melt model

xp   = 'IceOcean1rr';

for mi = 1:length(melt_models)
  melt = melt_models{mi};
  for ri = 1:length(resolutions)
    rres = ['r' resolutions{ri}];
    
    % Sliding laws
    MISOMIP1.(melt).(rres).slid = model;
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};

      if isempty(MISOMIP1.(xp).(rres).(slid).(melt).xGL); continue; end

      MISOMIP1.(melt).(rres).slid.nmod    =      MISOMIP1.(melt).(rres).slid.nmod + 1;
      MISOMIP1.(melt).(rres).slid.xGL_av  =      MISOMIP1.(melt).(rres).slid.xGL_av + MISOMIP1.(xp).(rres).(slid).(melt).xGL;
      MISOMIP1.(melt).(rres).slid.xGL_min = min( MISOMIP1.(melt).(rres).slid.xGL_min, MISOMIP1.(xp).(rres).(slid).(melt).xGL);
      MISOMIP1.(melt).(rres).slid.xGL_max = max( MISOMIP1.(melt).(rres).slid.xGL_max, MISOMIP1.(xp).(rres).(slid).(melt).xGL);

    end

    MISOMIP1.(melt).(rres).slid.xGL_av = MISOMIP1.(melt).(rres).slid.xGL_av / MISOMIP1.(melt).(rres).slid.nmod;

    % Sub-grid schemes
    MISOMIP1.(melt).(rres).sbg = model;
    for sbgi = 1:length(subgrid_schemes)
      sbg = subgrid_schemes{sbgi};

      if isempty(subgridmelt.(rres).(melt).(sbg).xGL); continue; end

      MISOMIP1.(melt).(rres).sbg.nmod    =      MISOMIP1.(melt).(rres).sbg.nmod + 1;
      MISOMIP1.(melt).(rres).sbg.xGL_av  =      MISOMIP1.(melt).(rres).sbg.xGL_av + subgridmelt.(rres).(melt).(sbg).xGL;
      MISOMIP1.(melt).(rres).sbg.xGL_min = min( MISOMIP1.(melt).(rres).sbg.xGL_min, subgridmelt.(rres).(melt).(sbg).xGL);
      MISOMIP1.(melt).(rres).sbg.xGL_max = max( MISOMIP1.(melt).(rres).sbg.xGL_max, subgridmelt.(rres).(melt).(sbg).xGL);

    end

    MISOMIP1.(melt).(rres).sbg.xGL_av = MISOMIP1.(melt).(rres).sbg.xGL_av / MISOMIP1.(melt).(rres).sbg.nmod;
  
  end
end

%% Plot plumes per experiment per melt model

wa = [400,170,400,170];
ha = [300];

margins_hor = [100,15,15,15,25];
margins_ver = [80,40];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + sum(wa);
hf = sum(margins_ver) + sum(ha);

H1.Fig = figure('position',[100,100,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    x = sum(margins_hor(1:i)) + sum(wa(1:i-1));
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*sum(ha(1:jp-1));
    H1.Ax(i,j) = axes('parent',H1.Fig,'units','pixels','position',[x,y,wa(i),ha(j)],'fontsize',24,...
      'xgrid','on','ygrid','on','xlim',[0,200]);
  end
end

H1.Ax = H1.Ax';

ylim = [300,475];
ytick = 300:25:475;

% Big axeses for plumes
title(H1.Ax(1,1),'A');
title(H1.Ax(1,3),'B');

set(H1.Ax(1,1),'ylim',ylim,'ytick',ytick);
set(H1.Ax(1,3),'ylim',ylim,'ytick',ytick,'yticklabels','');

set(H1.Ax(1,1),'xlim',[0,200],'xtick',0:50:200);
set(H1.Ax(1,3),'xlim',[0,200],'xtick',0:50:200);

ylabel(H1.Ax(1,1),'x_{GL} (km)')
xlabel(H1.Ax(1,1),'Time (yr)')
xlabel(H1.Ax(1,3),'Time (yr)')

% Small axeses for bars
set(H1.Ax(1,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H1.Ax(    1,2).XAxis.Visible = 'off';
H1.Ax(    1,2).YAxis.Visible = 'off';
set(H1.Ax(1,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H1.Ax(    1,4).XAxis.Visible = 'off';
H1.Ax(    1,4).YAxis.Visible = 'off';

% Empty line objects for legend
colors = lines(length(melt_models));
for mi = 1:length(melt_models)
  line('parent',H1.Ax(1,1),'xdata',[],'ydata',[],'linewidth',3,'color',colors(mi,:));
  line('parent',H1.Ax(1,2),'xdata',[],'ydata',[],'linewidth',3,'color',colors(mi,:));
end

% Separators in bar axeses
line('parent',H1.Ax(1,2),'xdata',[0,0]+0/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,2),'xdata',[0,0]+1/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,2),'xdata',[0,0]+2/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,2),'xdata',[0,0]+3/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,4),'xdata',[0,0]+0/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,4),'xdata',[0,0]+1/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,4),'xdata',[0,0]+2/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);
line('parent',H1.Ax(1,4),'xdata',[0,0]+3/3,'ydata',ylim,'color',[1,1,1]*0.7,'linewidth',1);

% Labels for bar axeses
dx = 0.1;
text(H1.Ax(1,2),0/3+dx,ylim(2)+range(ylim)*0.05,'0', 'fontsize',24);
text(H1.Ax(1,2),1/3+dx,ylim(2)+range(ylim)*0.05,'ra','fontsize',24);
text(H1.Ax(1,2),2/3+dx,ylim(2)+range(ylim)*0.05,'rr','fontsize',24);
text(H1.Ax(1,4),0/3+dx,ylim(2)+range(ylim)*0.05,'0', 'fontsize',24);
text(H1.Ax(1,4),1/3+dx,ylim(2)+range(ylim)*0.05,'ra','fontsize',24);
text(H1.Ax(1,4),2/3+dx,ylim(2)+range(ylim)*0.05,'rr','fontsize',24);

% Plot data
for mi = 1:length(melt_models)
  melt = melt_models{mi};
  color = colors(mi,:);

  % IceOcean0 (both axeses)
  plume = MISOMIP1.(melt).IceOcean0;
  xdata_lo =        plume.time;
  xdata_hi = flipud(plume.time);
  ydata_lo =        plume.xGL_min;
  ydata_hi = flipud(plume.xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  patch('parent',H1.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,1),'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  line('parent',H1.Ax(1,3),'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  
  % Bar
  nb   = 3 * (3*length(melt_models)+1);
  w    = 1/nb;
  x_lo = 0/3 + w + (mi-1)*3*w;
  x_hi = x_lo + 2*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  patch('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  line('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

  % IceOcean1r
  plume = MISOMIP1.(melt).IceOcean1ra;
  m = plume.time <= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,1),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);

  % IceOcean1ra
  plume = MISOMIP1.(melt).IceOcean1ra;
  m = plume.time >= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,1),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);
  
  % Bar
  nb   = 3 * (3*length(melt_models)+1);
  w    = 1/nb;
  x_lo = 1/3 + w + (mi-1)*3*w;
  x_hi = x_lo + 2*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

  % IceOcean1rr
  plume = MISOMIP1.(melt).IceOcean1rr;
  m = plume.time >= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,1),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);
  
  % Bar
  nb   = 3 * (3*length(melt_models)+1);
  w    = 1/nb;
  x_lo = 2/3 + w + (mi-1)*3*w;
  x_hi = x_lo + 2*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H1.Ax(1,2),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

  % IceOcean2r
  plume = MISOMIP1.(melt).IceOcean2ra;
  m = plume.time <= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,3),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);

  % IceOcean2ra
  plume = MISOMIP1.(melt).IceOcean2ra;
  m = plume.time >= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,3),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);
  
  % Bar
  nb   = 3 * (3*length(melt_models)+1);
  w    = 1/nb;
  x_lo = 1/3 + w + (mi-1)*3*w;
  x_hi = x_lo + 2*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

  % IceOcean2rr
  plume = MISOMIP1.(melt).IceOcean2rr;
  m = plume.time >= 100;
  xdata_lo =        plume.time(    m);
  xdata_hi = flipud(plume.time(    m));
  ydata_lo =        plume.xGL_min( m);
  ydata_hi = flipud(plume.xGL_max( m));
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H1.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H1.Ax(1,3),'xdata',plume.time( m),'ydata',plume.xGL_av( m)/1e3,'color',color,'linewidth',3);
  
  % Bar
  nb   = 3 * (3*length(melt_models)+1);
  w    = 1/nb;
  x_lo = 2/3 + w + (mi-1)*3*w;
  x_hi = x_lo + 2*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H1.Ax(1,4),'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  
end

% Legends
legend(H1.Ax(1,1),'Linear','Quadratic','M+','Plume','PICO','PICOP','location','southwest');

%% Plot sliding law/subgridscheme plumes per melt model

wa = [400,60,400,60,400,60];
ha = [300,300];

margins_hor = [100,15,15,15,15,15,25];
margins_ver = [80,40,40];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + sum(wa);
hf = sum(margins_ver) + sum(ha);

H2.Fig = figure('position',[200,200,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    x = sum(margins_hor(1:i)) + sum(wa(1:i-1));
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*sum(ha(1:jp-1));
    H2.Ax(i,j) = axes('parent',H2.Fig,'units','pixels','position',[x,y,wa(i),ha(j)],'fontsize',24,...
      'xgrid','on','ygrid','on','xlim',[0,200]);
  end
end

H2.Ax = H2.Ax';

ylim = [300,475];
ytick = 300:25:475;

% Big axeses for plumes
% title(H2.Ax(1,1),'A');
% title(H2.Ax(1,3),'B');
% title(H2.Ax(1,5),'C');
% title(H2.Ax(2,1),'D');
% title(H2.Ax(2,3),'E');
% title(H2.Ax(2,5),'F');
title(H2.Ax(1,1),'Linear');
title(H2.Ax(1,3),'Quadratic');
title(H2.Ax(1,5),'M+');
title(H2.Ax(2,1),'Plume');
title(H2.Ax(2,3),'PICO');
title(H2.Ax(2,5),'PICOP');

set(H2.Ax(1,1),'ylim',ylim,'ytick',ytick);
set(H2.Ax(1,3),'ylim',ylim,'ytick',ytick,'yticklabels','');
set(H2.Ax(1,5),'ylim',ylim,'ytick',ytick,'yticklabels','');
set(H2.Ax(2,1),'ylim',ylim,'ytick',ytick);
set(H2.Ax(2,3),'ylim',ylim,'ytick',ytick,'yticklabels','');
set(H2.Ax(2,5),'ylim',ylim,'ytick',ytick,'yticklabels','');

set(H2.Ax(1,1),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H2.Ax(1,3),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H2.Ax(1,5),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H2.Ax(2,1),'xlim',[0,200],'xtick',0:50:200);
set(H2.Ax(2,3),'xlim',[0,200],'xtick',0:50:200);
set(H2.Ax(2,5),'xlim',[0,200],'xtick',0:50:200);

ylabel(H2.Ax(1,1),'x_{GL} (km)')
ylabel(H2.Ax(2,1),'x_{GL} (km)')
xlabel(H2.Ax(2,1),'Time (yr)')
xlabel(H2.Ax(2,3),'Time (yr)')
xlabel(H2.Ax(2,5),'Time (yr)')

% Small axeses for bars
set(H2.Ax(1,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    1,2).XAxis.Visible = 'off';
H2.Ax(    1,2).YAxis.Visible = 'off';
set(H2.Ax(1,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    1,4).XAxis.Visible = 'off';
H2.Ax(    1,4).YAxis.Visible = 'off';
set(H2.Ax(1,6),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    1,6).XAxis.Visible = 'off';
H2.Ax(    1,6).YAxis.Visible = 'off';
set(H2.Ax(2,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    2,2).XAxis.Visible = 'off';
H2.Ax(    2,2).YAxis.Visible = 'off';
set(H2.Ax(2,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    2,4).XAxis.Visible = 'off';
H2.Ax(    2,4).YAxis.Visible = 'off';
set(H2.Ax(2,6),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'ygrid','on','xgrid','off');
H2.Ax(    2,6).XAxis.Visible = 'off';
H2.Ax(    2,6).YAxis.Visible = 'off';

% Empty objects for legends
colors = lines(4);
patch('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(1,:),'facealpha',0.45,'edgecolor','none');
patch('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(2,:),'facealpha',0.45,'edgecolor','none');
patch('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(3,:),'facealpha',0.45,'edgecolor','none');
patch('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(4,:),'facealpha',0.45,'edgecolor','none');

% Plot data
axlin  = [H2.Ax(1,1), H2.Ax(1,3), H2.Ax(1,5), H2.Ax(2,1), H2.Ax(2,3), H2.Ax(2,5)];
axlin2 = [H2.Ax(1,2), H2.Ax(1,4), H2.Ax(1,6), H2.Ax(2,2), H2.Ax(2,4), H2.Ax(2,6)];
for mi = 1:length(melt_models)
  melt = melt_models{mi};
  ax   = axlin(mi);
  ax2  = axlin2(mi);
  
  % Sliding law plumes - 5 km
  color = colors(1,:);
  plume = MISOMIP1.(melt).r5km.slid;
  xdata_lo =        plume.time;
  xdata_hi = flipud(plume.time);
  ydata_lo =        plume.xGL_min;
  ydata_hi = flipud(plume.xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',ax,'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',ax,'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  
  % Sliding law bars - 5 km
  nb   = 13;
  w    = 1/nb;
  x_lo = 1*w;
  x_hi = 3*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',ax2,'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',ax2,'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  
  % Sliding law plumes - 2 km
  color = colors(2,:);
  plume = MISOMIP1.(melt).r2km.slid;
  xdata_lo =        plume.time;
  xdata_hi = flipud(plume.time);
  ydata_lo =        plume.xGL_min;
  ydata_hi = flipud(plume.xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',ax,'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',ax,'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  
  % Sliding law bars - 2 km
  nb   = 13;
  w    = 1/nb;
  x_lo = 4*w;
  x_hi = 6*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',ax2,'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',ax2,'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  
  % Sub-grid scheme plumes - 5 km
  color = colors(3,:);
  plume = MISOMIP1.(melt).r5km.sbg;
  xdata_lo =        plume.time;
  xdata_hi = flipud(plume.time);
  ydata_lo =        plume.xGL_min;
  ydata_hi = flipud(plume.xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',ax,'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',ax,'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  
  % Sub-grid scheme bars - 5 km
  nb   = 13;
  w    = 1/nb;
  x_lo = 7*w;
  x_hi = 9*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',ax2,'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',ax2,'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  
  % Sub-grid scheme plumes - 2 km
  color = colors(4,:);
  plume = MISOMIP1.(melt).r2km.sbg;
  xdata_lo =        plume.time;
  xdata_hi = flipud(plume.time);
  ydata_lo =        plume.xGL_min;
  ydata_hi = flipud(plume.xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',ax,'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',ax,'xdata',plume.time,'ydata',plume.xGL_av/1e3,'color',color,'linewidth',3);
  
  % Sub-grid scheme bars - 2 km
  nb   = 13;
  w    = 1/nb;
  x_lo = 10*w;
  x_hi = 12*w;
  y_lo = plume.xGL_min( end);
  y_av = plume.xGL_av(  end);
  y_hi = plume.xGL_max( end);
  patch('parent',ax2,'xdata',[x_lo,x_hi,x_hi,x_lo],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,...
    'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',ax2,'xdata',[x_lo,x_hi],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
  
end

% Legend
legend(H2.Ax(1,1),'Sliding law - 5km','Sliding law - 2km','Sub-grid scheme - 5 km','Sub-grid scheme - 2 km','location','southwest');