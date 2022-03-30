clc
clear all
close all

%% Collect data from Cornford et al. 2020

if exist('Cornford_etal_2020_data.mat','file')
  load('Cornford_etal_2020_data.mat')
else

  filenames = dir('/Users/berends/Documents/Datasets/MISMIP+/supplement/submission_data');
  
  Cornford.ice1r.models  = [];
  Cornford.ice1ra.models = [];
  Cornford.ice1rr.models = [];
  Cornford.ice2r.models  = [];
  Cornford.ice2ra.models = [];
  Cornford.ice2rr.models = [];
  
  % Read data
  for fi = 1:length(filenames)
    filename = filenames(fi).name;
    if length(filename)<3; continue; end
    if strcmpi(filename(end-2:end),'.nc')

      if strcmpi(filename(1:4),'Ice0')
        exp = 'ice0';
        model.Name = filename(6:end-3);
        % Don't do anything with ice0...
        continue
      elseif strcmpi(filename(1:6),'Ice1ra')
        exp = 'ice1ra';
        model.Name = filename(8:end-3);
      elseif strcmpi(filename(1:6),'Ice1rr')
        exp = 'ice1rr';
        model.Name = filename(8:end-3);
      elseif strcmpi(filename(1:6),'Ice2ra')
        exp = 'ice2ra';
        model.Name = filename(8:end-3);
      elseif strcmpi(filename(1:6),'Ice2rr')
        exp = 'ice2rr';
        model.Name = filename(8:end-3);
      elseif strcmpi(filename(1:5),'Ice1r')
        exp = 'ice1r';
        model.Name = filename(7:end-3);
      elseif strcmpi(filename(1:5),'Ice2r')
        exp = 'ice2r';
        model.Name = filename(7:end-3);
      else
        error('unrecognised MISMIP+ experiment')
      end
      
      disp(['Reading data from Cornford et al. (2020), experiment ' exp ', model "' model.Name '"'])

      filename_full = ['/Users/berends/Documents/Datasets/MISMIP+/supplement/submission_data/' filename];

      % Read data, get xGL at mid-channel
      time = ncread(filename_full,'time')';
      xGLf = ncread(filename_full,'xGL');
      yGLf = ncread(filename_full,'yGL');
      
      % Mark NaN values
      yGLf( yGLf>8e5) = NaN;
      xGLf( xGLf>8e5) = NaN;
      
      % Get everything in [m]
      if max(xGLf(:)<800)
        xGLf = xGLf * 1e3;
        yGLf = yGLf * 1e3;
      end
      
      % Get mid-channel xGL
      xGL = zeros(size(time));
      for ti = 1:length(time)
        yGLti = yGLf( ti,:);
        dy = abs(yGLti - 40000);
        j = find(dy==min(dy));
        if isempty(j); xGL(ti)=NaN; continue; end
        j=j(1);
        xGL(ti) = xGLf(ti,j);
      end
      
      % Remove values where time is NaN (why even have this?)
      xGL( isnan(time)) = [];
      time(isnan(time)) = [];
      
      % If data is bullshit, skip it
%       if sum(isnan(xGL))>0; continue; end
      hasdoubletimes = false;
      for ti1 = 1:length(time)
        for ti2 = ti1+1:length(time)
          if time(ti2)==time(ti1)
            hasdoubletimes = true;
          end
        end
      end
      if hasdoubletimes; continue; end
      
      % Interpolate to proper grid
      if strcmpi(exp,'ice1r') || strcmpi(exp,'ice2r')
        model.time = (0:10:100)';
      else
        model.time = (100:10:200)';
      end
      model.xGL = interp1(time,xGL,model.time);
      
      Cornford.(exp).models{end+1} = model;

    end
  end
  
  %% Only keep the "main subset", i.e. those models that:
  %   1) have completed all experiments
  %   2) place their Ice0 steady-state grounding line at x = 450 ± 15 km2
  %   3) see their grounding lines retreat along the center of the channel to x = 385 ± 35 km.
  
  experiments = {'ice1r','ice1ra','ice1rr','ice2r','ice2ra','ice2rr'};
  
  % List all participating models
  models = {};
  for xi = 1:length(experiments)
    exp = experiments{xi};
    for mi = 1:length(Cornford.(exp).models)
      is_listed = false;
      for mii = 1:length(models)
        if strcmpi(models{mii}.Name,Cornford.(exp).models{mi}.Name)
          is_listed = true;
          break
        end
      end
      if ~is_listed
        mii = length(models)+1;
        models{mii}.Name      = Cornford.(exp).models{mi}.Name;
        models{mii}.ice1r     = false;
        models{mii}.ice1ra    = false;
        models{mii}.ice1rr    = false;
        models{mii}.ice2r     = false;
        models{mii}.ice2ra    = false;
        models{mii}.ice2rr    = false;
        models{mii}.xGL_ice0  = Inf;
        models{mii}.xGL_ice1r = Inf;
        models{mii}.isvalid   = true;
      end
    end
  end
  
  % For each model, list which experiments they have completed
  for xi = 1:length(experiments)
    exp = experiments{xi};
    for mi = 1:length(Cornford.(exp).models)
      for mii = 1:length(models)
        if strcmpi(models{mii}.Name,Cornford.(exp).models{mi}.Name)
          models{mii}.(exp) = true;
        end
      end
    end
  end
  
  % List GL values for all models
  for mii = 1:length(models)
    for mi = 1:length(Cornford.ice1r.models)
      if strcmpi(Cornford.ice1r.models{mi}.Name,models{mii}.Name)
        models{mii}.xGL_ice0  = Cornford.ice1r.models{mi}.xGL(1);
        models{mii}.xGL_ice1r = Cornford.ice1r.models{mi}.xGL(end);
      end
    end
  end
  
  % Mark invalid models
  for mii = 1:length(models)
    if ~models{mii}.ice1r || ~models{mii}.ice1ra || ~models{mii}.ice1rr || ...
       ~models{mii}.ice2r || ~models{mii}.ice2ra || ~models{mii}.ice2rr || ...
       abs(models{mii}.xGL_ice0/1e3-450)>15 || abs(models{mii}.xGL_ice1r/1e3-385)>35
     models{mii}.isvalid = false;
    end
  end
  
  % Remove invalid models
  for xi = 1:length(experiments)
    exp = experiments{xi};
    mi = 1;
    while mi<=length(Cornford.(exp).models)
      
      % Check if this model is valid
      isvalid = false;
      for mii = 1:length(models)
        if strcmpi(models{mii}.Name,Cornford.(exp).models{mi}.Name)
          isvalid = models{mii}.isvalid;
        end
      end
      if ~isvalid
        Cornford.(exp).models(mi) = [];
      else
        mi = mi+1;
      end
    end
  end
  
  % Calculate plumes
  Cornford.ice1r.time  =   0:10:100;
  Cornford.ice1ra.time = 100:10:200;
  Cornford.ice1rr.time = 100:10:200;
  Cornford.ice2r.time  =   0:10:100;
  Cornford.ice2ra.time = 100:10:200;
  Cornford.ice2rr.time = 100:10:200;
  
  for xi = 1:length(experiments)
    
    exp = experiments{xi};
    
    nmodels = zeros(1,11);
    
    Cornford.(exp).xGL_av  = zeros(size(Cornford.(exp).time));
    Cornford.(exp).xGL_min = zeros(size(Cornford.(exp).time)) + Inf;
    Cornford.(exp).xGL_max = zeros(size(Cornford.(exp).time)) - Inf;
    
    for mi = 1:length(Cornford.(exp).models)
      
      for ti = 1:11
        x = Cornford.(exp).models{mi}.xGL(ti);
        if ~isnan(x) && x > 0 && x < 8e5
          nmodels( ti) = nmodels( ti) + 1;
          Cornford.(exp).xGL_av(  ti)  =     Cornford.(exp).xGL_av( ti) + x;
          Cornford.(exp).xGL_min( ti) = min(Cornford.(exp).xGL_min( ti), x);
          Cornford.(exp).xGL_max( ti) = max(Cornford.(exp).xGL_max( ti), x);
        end
      
      end
    end
    
    Cornford.(exp).xGL_av = Cornford.(exp).xGL_av ./ nmodels;
    
  end
  
  save('Cornford_etal_2020_data.mat','Cornford')
  
end % if exist('Cornford_etal_2020_data.mat','file')

%% Collect data from IMAU-ICE

experiments     = {'ice0','ice1r','ice1ra','ice1rr','ice2r','ice2ra','ice2rr'};
resolutions     = {'2km','5km'};
sliding_laws    = {'Weertman','Tsai2015','Schoof2005'};
stress_balances = {'DIVA','hybrid'};
subgrid_schemes = {'FCMP','PMP','NMP'};

IMAU_ICE.nruns   = 0;
IMAU_ICE.results = [];

for i = 1:length(experiments    ); IMAU_ICE.(['ens_' experiments{     i}]).list = []; end
for i = 1:length(resolutions    ); IMAU_ICE.(['ens_' resolutions{     i}]).list = []; end
for i = 1:length(sliding_laws   ); IMAU_ICE.(['ens_' sliding_laws{    i}]).list = []; end
for i = 1:length(stress_balances); IMAU_ICE.(['ens_' stress_balances{ i}]).list = []; end
for i = 1:length(subgrid_schemes); IMAU_ICE.(['ens_' subgrid_schemes{ i}]).list = []; end

henk = dir;
for i = 1:length(henk)
  
  if henk(i).isdir && contains(henk(i).name,'MISMIPplus_ice')
    
    IMAU_ICE.nruns = IMAU_ICE.nruns + 1;
    
    filename = [henk(i).name '/MISMIPplus_output.nc'];
    
    time = ncread(filename,'time');
    xGL  = ncread(filename,'xGL');
    ny  = size(xGL,2);
    jmid = ceil(ny/2);
    xGL = xGL(:,jmid);
    
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
    IMAU_ICE.(['ens_' experiment    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' resolution    ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' sliding_law   ]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' stress_balance]).list(end+1) = IMAU_ICE.nruns;
    IMAU_ICE.(['ens_' subgrid_scheme]).list(end+1) = IMAU_ICE.nruns;
    
  end
  
end

%% Create plumes per experiment
for xi = 1:length(experiments)
  
  list    = IMAU_ICE.(['ens_' experiments{xi}]).list;
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
  
  IMAU_ICE.(['ens_' experiments{xi}]).time    = time;
  IMAU_ICE.(['ens_' experiments{xi}]).xGL_av  = xGL_av;
  IMAU_ICE.(['ens_' experiments{xi}]).xGL_min = xGL_min;
  IMAU_ICE.(['ens_' experiments{xi}]).xGL_max = xGL_max;
end

%% Compare to results from Cornford et al. (2020)

wa = 400;
ha = 300;

margins_hor = [100,45,90];
margins_ver = [80,40];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H1.Fig = figure('position',[100,100,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    x = sum(margins_hor(1:i)) + (i-1)*wa;
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*ha;
    H1.Ax(i,j) = axes('parent',H1.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24,...
      'xgrid','on','ygrid','on','xlim',[0,200]);
  end
end

H1.Ax = H1.Ax';

set(H1.Ax(1,1),'ylim',[300,475],'ytick',300:25:475);
set(H1.Ax(1,2),'ylim',[340,460],'ytick',340:20:460,'yaxislocation','right');

xlabel(H1.Ax(1,1),'Time (yr)')
xlabel(H1.Ax(1,2),'Time (yr)')
ylabel(H1.Ax(1,1),'x_{GL} (km)')
ylabel(H1.Ax(1,2),'x_{GL} (km)')

title(H1.Ax(1,1),'A')
title(H1.Ax(1,2),'B')

% Empty line objects for legend
color = [0,0,0];
patch('parent',H1.Ax(1,1),'xdata',[],'ydata',[],'facecolor',color,'facealpha',0.45,'edgecolor','none');
color = [0,0,1];
patch('parent',H1.Ax(1,1),'xdata',[],'ydata',[],'facecolor',color,'facealpha',0.45,'edgecolor','none');
color = [1,0,0];
patch('parent',H1.Ax(1,1),'xdata',[],'ydata',[],'facecolor',color,'facealpha',0.45,'edgecolor','none');
color = [1.0,0.63,0.0];
patch('parent',H1.Ax(1,1),'xdata',[],'ydata',[],'facecolor',color,'facealpha',0.45,'edgecolor','none');

line('parent',H1.Ax(1,2),'xdata',[],'ydata',[],'linewidth',3,'color','k');
patch('parent',H1.Ax(1,2),'xdata',[],'ydata',[],'facecolor',[0,0,0],'facealpha',0.45,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',[],'ydata',[],'linewidth',2,'color','k','linestyle',':');
line('parent',H1.Ax(1,2),'xdata',[],'ydata',[],'linewidth',2,'color','k','linestyle','--');

% Plot Cornford et al. (2020) data with dashed lines showing min-max ranges

% ice1r
color = [0,0,1];
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1r.time,'ydata',Cornford.ice1r.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1r.time,'ydata',Cornford.ice1r.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1r.time,'ydata',Cornford.ice1r.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');
% ice1ra
color = [1,0,0];
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1ra.time,'ydata',Cornford.ice1ra.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1ra.time,'ydata',Cornford.ice1ra.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1ra.time,'ydata',Cornford.ice1ra.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');
% ice1rr
color = [1.0,0.63,0.0];
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1rr.time,'ydata',Cornford.ice1rr.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1rr.time,'ydata',Cornford.ice1rr.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,1),'xdata',Cornford.ice1rr.time,'ydata',Cornford.ice1rr.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');
% ice2r
color = [0,0,1];
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2r.time,'ydata',Cornford.ice2r.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2r.time,'ydata',Cornford.ice2r.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2r.time,'ydata',Cornford.ice2r.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');
% ice2ra
color = [1,0,0];
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2ra.time,'ydata',Cornford.ice2ra.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2ra.time,'ydata',Cornford.ice2ra.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2ra.time,'ydata',Cornford.ice2ra.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');
% ice2rr
color = [1.0,0.63,0.0];
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2rr.time,'ydata',Cornford.ice2rr.xGL_min/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2rr.time,'ydata',Cornford.ice2rr.xGL_max/1e3,'color',color,'linewidth',2,'linestyle','--');
line('parent',H1.Ax(1,2),'xdata',Cornford.ice2rr.time,'ydata',Cornford.ice2rr.xGL_av /1e3,'color',color,'linewidth',2,'linestyle',':');

% Plot IMAU-ICE simulations as plumes

% ice0 (both axeses)
color = [0.0,0.0,0.0];
xdata_lo =        IMAU_ICE.ens_ice0.time;
xdata_hi = flipud(IMAU_ICE.ens_ice0.time);
ydata_lo =        IMAU_ICE.ens_ice0.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice0.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ens_ice0.time,'ydata',IMAU_ICE.ens_ice0.xGL_av/1e3,'color',color,'linewidth',3);
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ens_ice0.time,'ydata',IMAU_ICE.ens_ice0.xGL_av/1e3,'color',color,'linewidth',3);

% ice1r
color = [0,0,1];
xdata_lo =        IMAU_ICE.ens_ice1r.time;
xdata_hi = flipud(IMAU_ICE.ens_ice1r.time);
ydata_lo =        IMAU_ICE.ens_ice1r.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice1r.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ens_ice1r.time,'ydata',IMAU_ICE.ens_ice1r.xGL_av/1e3,'color',color,'linewidth',3);

% ice1ra
color = [1,0,0];
xdata_lo =        IMAU_ICE.ens_ice1ra.time;
xdata_hi = flipud(IMAU_ICE.ens_ice1ra.time);
ydata_lo =        IMAU_ICE.ens_ice1ra.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice1ra.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ens_ice1ra.time,'ydata',IMAU_ICE.ens_ice1ra.xGL_av/1e3,'color',color,'linewidth',3);

% ice1rr
color = [1.0,0.63,0.0];
xdata_lo =        IMAU_ICE.ens_ice1rr.time;
xdata_hi = flipud(IMAU_ICE.ens_ice1rr.time);
ydata_lo =        IMAU_ICE.ens_ice1rr.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice1rr.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ens_ice1rr.time,'ydata',IMAU_ICE.ens_ice1rr.xGL_av/1e3,'color',color,'linewidth',3);

% ice2r
color = [0,0,1];
xdata_lo =        IMAU_ICE.ens_ice2r.time;
xdata_hi = flipud(IMAU_ICE.ens_ice2r.time);
ydata_lo =        IMAU_ICE.ens_ice2r.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice2r.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ens_ice2r.time,'ydata',IMAU_ICE.ens_ice2r.xGL_av/1e3,'color',color,'linewidth',3);

% ice2ra
color = [1,0,0];
m = IMAU_ICE.ens_ice2ra.time >= 100;
xdata_lo =        IMAU_ICE.ens_ice2ra.time(    m);
xdata_hi = flipud(IMAU_ICE.ens_ice2ra.time(    m));
ydata_lo =        IMAU_ICE.ens_ice2ra.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ens_ice2ra.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ens_ice2ra.time,'ydata',IMAU_ICE.ens_ice2ra.xGL_av/1e3,'color',color,'linewidth',3);

% ice2rr
color = [1.0,0.63,0.0];
xdata_lo =        IMAU_ICE.ens_ice2rr.time;
xdata_hi = flipud(IMAU_ICE.ens_ice2rr.time);
ydata_lo =        IMAU_ICE.ens_ice2rr.xGL_min;
ydata_hi = flipud(IMAU_ICE.ens_ice2rr.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ens_ice2rr.time,'ydata',IMAU_ICE.ens_ice2rr.xGL_av/1e3,'color',color,'linewidth',3);

% Legends
legend(H1.Ax(1,1),'ice0','ice1r','ice1ra','ice1rr','location','southwest');
legend(H1.Ax(1,2),'IMAU-ICE mean','IMAU-ICE range','C20 mean','C20 range','location','southwest');
