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

experiments     = {'ice0','ice1ra','ice1rr','ice2ra','ice2rr'};
resolutions     = {'5km','2km'};
sliding_laws    = {'slid1','slid2','slid3'};
stress_balances = {'DIVA','SIASSA'};
subgrid_schemes = {'NMP','FCMP','PMP'};

for xi = 1:length(experiments)
  xp = experiments{xi};
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};
      for sti = 1:length(stress_balances)
        strb = stress_balances{sti};
        for sbgi = 1:length(subgrid_schemes)
          sbg = subgrid_schemes{sbgi};
      
          foldername = ['MISMIPplus_' xp '_' res '_' slid '_' strb '_' sbg];
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

          IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time = time;
          IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL  = xGL;
      
        end
      end
    end
  end
end

%% Create plumes per experiment

for xi = 1:length(experiments)
  xp = experiments{xi};
  
  time    = (0:5:200)';
  xGL_av  = zeros(size(time));
  xGL_min = zeros(size(time)) + Inf;
  xGL_max = zeros(size(time)) - Inf;
  nmod    = 0;
  
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};
      for sti = 1:length(stress_balances)
        strb = stress_balances{sti};
        for sbgi = 1:length(subgrid_schemes)
          sbg = subgrid_schemes{sbgi};
          
          if ~isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time)
            nmod = nmod + 1;
            xGL_av  =      xGL_av + IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
            xGL_min = min( xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
            xGL_max = max( xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          end
          
        end
      end
    end
  end
  
  if nmod > 0
    xGL_av = xGL_av / nmod;
  end
  
  IMAU_ICE.(xp).time    = time;
  IMAU_ICE.(xp).xGL_av  = xGL_av;
  IMAU_ICE.(xp).xGL_min = xGL_min;
  IMAU_ICE.(xp).xGL_max = xGL_max;
  
end

%% Create plumes per resolution for ice1rr

xp = 'ice1rr';

for ri = 1:length(resolutions)
  res  = resolutions{ri};
  rres = ['r' resolutions{ri}];
  
  time    = (0:5:200)';
  xGL_av  = zeros(size(time));
  xGL_min = zeros(size(time)) + Inf;
  xGL_max = zeros(size(time)) - Inf;
  nmod    = 0;
  
  for sli = 1:length(sliding_laws)
    slid = sliding_laws{sli};
    for sti = 1:length(stress_balances)
      strb = stress_balances{sti};
      for sbgi = 1:length(subgrid_schemes)
        sbg = subgrid_schemes{sbgi};

        if ~isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time)
          nmod = nmod + 1;
          xGL_av  =      xGL_av + IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          xGL_min = min( xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          xGL_max = max( xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
        end

      end
    end
  end
  
  if nmod > 0
    xGL_av = xGL_av / nmod;
  end
  
  IMAU_ICE.(rres).time    = time;
  IMAU_ICE.(rres).xGL_av  = xGL_av;
  IMAU_ICE.(rres).xGL_min = xGL_min;
  IMAU_ICE.(rres).xGL_max = xGL_max;
  
end

%% Create plumes per sliding law for ice1rr

xp = 'ice1rr';

for sli = 1:length(sliding_laws)
  slid = sliding_laws{sli};
  
  time    = (0:5:200)';
  xGL_av  = zeros(size(time));
  xGL_min = zeros(size(time)) + Inf;
  xGL_max = zeros(size(time)) - Inf;
  nmod    = 0;
  
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sti = 1:length(stress_balances)
      strb = stress_balances{sti};
      for sbgi = 1:length(subgrid_schemes)
        sbg = subgrid_schemes{sbgi};

        if ~isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time)
          nmod = nmod + 1;
          xGL_av  =      xGL_av + IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          xGL_min = min( xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          xGL_max = max( xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
        end

      end
    end
  end
  
  if nmod > 0
    xGL_av = xGL_av / nmod;
  end
  
  IMAU_ICE.(slid).time    = time;
  IMAU_ICE.(slid).xGL_av  = xGL_av;
  IMAU_ICE.(slid).xGL_min = xGL_min;
  IMAU_ICE.(slid).xGL_max = xGL_max;
  
end

%% Create plumes per stress balance for ice1rr

xp = 'ice1rr';

for sti = 1:length(stress_balances)
  strb = stress_balances{sti};
  
  time    = (0:5:200)';
  xGL_av  = zeros(size(time));
  xGL_min = zeros(size(time)) + Inf;
  xGL_max = zeros(size(time)) - Inf;
  nmod    = 0;
  
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};
      for sbgi = 1:length(subgrid_schemes)
        sbg = subgrid_schemes{sbgi};

        if ~isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time)
          nmod = nmod + 1;
          xGL_av  =      xGL_av + IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          xGL_min = min( xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          xGL_max = max( xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
        end

      end
    end
  end
  
  if nmod > 0
    xGL_av = xGL_av / nmod;
  end
  
  IMAU_ICE.(strb).time    = time;
  IMAU_ICE.(strb).xGL_av  = xGL_av;
  IMAU_ICE.(strb).xGL_min = xGL_min;
  IMAU_ICE.(strb).xGL_max = xGL_max;
  
end

%% Create plumes per subgrid scheme for ice1rr

xp = 'ice1rr';

for sbgi = 1:length(subgrid_schemes)
  sbg = subgrid_schemes{sbgi};
  
  time    = (0:5:200)';
  xGL_av  = zeros(size(time));
  xGL_min = zeros(size(time)) + Inf;
  xGL_max = zeros(size(time)) - Inf;
  nmod    = 0;
  
  for ri = 1:length(resolutions)
    res  = resolutions{ri};
    rres = ['r' resolutions{ri}];
    for sli = 1:length(sliding_laws)
      slid = sliding_laws{sli};
      for sti = 1:length(stress_balances)
        strb = stress_balances{sti};

        if ~isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).time)
          nmod = nmod + 1;
          xGL_av  =      xGL_av + IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          xGL_min = min( xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          xGL_max = max( xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
        end

      end
    end
  end
  
  if nmod > 0
    xGL_av = xGL_av / nmod;
  end
  
  IMAU_ICE.(sbg).time    = time;
  IMAU_ICE.(sbg).xGL_av  = xGL_av;
  IMAU_ICE.(sbg).xGL_min = xGL_min;
  IMAU_ICE.(sbg).xGL_max = xGL_max;
  
end

%% Create plumes for benchmark run
xp = 'ice1rr';
benchmark_5km.rres.name = 'r5km';
benchmark_5km.slid.name = 'slid2';
benchmark_5km.strb.name = 'DIVA';
benchmark_5km.sbg.name  = 'FCMP';

time = (0:5:200)';

benchmark_5km.rres.time    = time;
benchmark_5km.rres.xGL     = zeros(size(time));
benchmark_5km.rres.xGL_min = zeros(size(time)) + Inf;
benchmark_5km.rres.xGL_max = zeros(size(time)) - Inf;

benchmark_5km.slid.time    = time;
benchmark_5km.slid.xGL     = zeros(size(time));
benchmark_5km.slid.xGL_min = zeros(size(time)) + Inf;
benchmark_5km.slid.xGL_max = zeros(size(time)) - Inf;

benchmark_5km.strb.time    = time;
benchmark_5km.strb.xGL     = zeros(size(time));
benchmark_5km.strb.xGL_min = zeros(size(time)) + Inf;
benchmark_5km.strb.xGL_max = zeros(size(time)) - Inf;

benchmark_5km.sbg.time     = time;
benchmark_5km.sbg.xGL      = zeros(size(time));
benchmark_5km.sbg.xGL_min  = zeros(size(time)) + Inf;
benchmark_5km.sbg.xGL_max  = zeros(size(time)) - Inf;

benchmark_2km = benchmark_5km;
benchmark_2km.rres.name = 'r2km';

for ri = 1:length(resolutions)
  res  = resolutions{ri};
  rres = ['r' resolutions{ri}];
  for sli = 1:length(sliding_laws)
    slid = sliding_laws{sli};
    for sti = 1:length(stress_balances)
      strb = stress_balances{sti};
      for sbgi = 1:length(subgrid_schemes)
        sbg = subgrid_schemes{sbgi};
        
        if isempty(IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL); continue; end
        
        %% 5 km baseline
        
        % Resolution
        if  strcmpi(slid,benchmark_5km.slid.name) && ...
            strcmpi(strb,benchmark_5km.strb.name) && ...
            strcmpi(sbg, benchmark_5km.sbg.name)
          benchmark_5km.rres.xGL_min = min( benchmark_5km.rres.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_5km.rres.xGL_max = max( benchmark_5km.rres.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(rres,benchmark_5km.rres.name)
            benchmark_5km.rres.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Sliding law
        if  strcmpi(rres,benchmark_5km.rres.name) && ...
            strcmpi(strb,benchmark_5km.strb.name) && ...
            strcmpi(sbg, benchmark_5km.sbg.name)
          benchmark_5km.slid.xGL_min = min( benchmark_5km.slid.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_5km.slid.xGL_max = max( benchmark_5km.slid.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(slid,benchmark_5km.slid.name)
            benchmark_5km.slid.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Stress balance
        if  strcmpi(rres,benchmark_5km.rres.name) && ...
            strcmpi(slid,benchmark_5km.slid.name) && ...
            strcmpi(sbg, benchmark_5km.sbg.name)
          benchmark_5km.strb.xGL_min = min( benchmark_5km.strb.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_5km.strb.xGL_max = max( benchmark_5km.strb.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(strb,benchmark_5km.strb.name)
            benchmark_5km.strb.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Sub-grid scheme
        if  strcmpi(rres,benchmark_5km.rres.name) && ...
            strcmpi(slid,benchmark_5km.slid.name) && ...
            strcmpi(strb,benchmark_5km.strb.name)
          benchmark_5km.sbg.xGL_min  = min( benchmark_5km.sbg.xGL_min,  IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_5km.sbg.xGL_max  = max( benchmark_5km.sbg.xGL_max,  IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(sbg,benchmark_5km.sbg.name)
            benchmark_5km.sbg.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        %% 2 km baseline
        
        % Resolution
        if  strcmpi(slid,benchmark_2km.slid.name) && ...
            strcmpi(strb,benchmark_2km.strb.name) && ...
            strcmpi(sbg, benchmark_2km.sbg.name)
          benchmark_2km.rres.xGL_min = min( benchmark_2km.rres.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_2km.rres.xGL_max = max( benchmark_2km.rres.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(rres,benchmark_2km.rres.name)
            benchmark_2km.rres.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Sliding law
        if  strcmpi(rres,benchmark_2km.rres.name) && ...
            strcmpi(strb,benchmark_2km.strb.name) && ...
            strcmpi(sbg, benchmark_2km.sbg.name)
          benchmark_2km.slid.xGL_min = min( benchmark_2km.slid.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_2km.slid.xGL_max = max( benchmark_2km.slid.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(slid,benchmark_2km.slid.name)
            benchmark_2km.slid.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Stress balance
        if  strcmpi(rres,benchmark_2km.rres.name) && ...
            strcmpi(slid,benchmark_2km.slid.name) && ...
            strcmpi(sbg, benchmark_2km.sbg.name)
          benchmark_2km.strb.xGL_min = min( benchmark_2km.strb.xGL_min, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_2km.strb.xGL_max = max( benchmark_2km.strb.xGL_max, IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(strb,benchmark_2km.strb.name)
            benchmark_2km.strb.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
        % Sub-grid scheme
        if  strcmpi(rres,benchmark_2km.rres.name) && ...
            strcmpi(slid,benchmark_2km.slid.name) && ...
            strcmpi(strb,benchmark_2km.strb.name)
          benchmark_2km.sbg.xGL_min  = min( benchmark_2km.sbg.xGL_min,  IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          benchmark_2km.sbg.xGL_max  = max( benchmark_2km.sbg.xGL_max,  IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL);
          if strcmpi(sbg,benchmark_2km.sbg.name)
            benchmark_2km.sbg.xGL = IMAU_ICE.(xp).(rres).(slid).(strb).(sbg).xGL;
          end
        end
        
      end
    end
  end
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
xdata_lo =        IMAU_ICE.ice0.time;
xdata_hi = flipud(IMAU_ICE.ice0.time);
ydata_lo =        IMAU_ICE.ice0.xGL_min;
ydata_hi = flipud(IMAU_ICE.ice0.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ice0.time,'ydata',IMAU_ICE.ice0.xGL_av/1e3,'color',color,'linewidth',3);
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ice0.time,'ydata',IMAU_ICE.ice0.xGL_av/1e3,'color',color,'linewidth',3);

% ice1r
color = [0,0,1];
m = IMAU_ICE.ice1ra.time <= 100;
xdata_lo =        IMAU_ICE.ice1ra.time(    m);
xdata_hi = flipud(IMAU_ICE.ice1ra.time(    m));
ydata_lo =        IMAU_ICE.ice1ra.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice1ra.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ice1ra.time( m),'ydata',IMAU_ICE.ice1ra.xGL_av( m)/1e3,'color',color,'linewidth',3);

% ice1ra
color = [1,0,0];
m = IMAU_ICE.ice1ra.time >= 100;
xdata_lo =        IMAU_ICE.ice1ra.time(    m);
xdata_hi = flipud(IMAU_ICE.ice1ra.time(    m));
ydata_lo =        IMAU_ICE.ice1ra.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice1ra.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ice1ra.time( m),'ydata',IMAU_ICE.ice1ra.xGL_av( m)/1e3,'color',color,'linewidth',3);

% ice1rr
color = [1.0,0.63,0.0];
m = IMAU_ICE.ice1rr.time >= 100;
xdata_lo =        IMAU_ICE.ice1rr.time(    m);
xdata_hi = flipud(IMAU_ICE.ice1rr.time(    m));
ydata_lo =        IMAU_ICE.ice1rr.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice1rr.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,1),'xdata',IMAU_ICE.ice1rr.time( m),'ydata',IMAU_ICE.ice1rr.xGL_av( m)/1e3,'color',color,'linewidth',3);

% ice2r
color = [0,0,1];
m = IMAU_ICE.ice2ra.time <= 100;
xdata_lo =        IMAU_ICE.ice2ra.time(    m);
xdata_hi = flipud(IMAU_ICE.ice2ra.time(    m));
ydata_lo =        IMAU_ICE.ice2ra.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice2ra.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ice2ra.time( m),'ydata',IMAU_ICE.ice2ra.xGL_av( m)/1e3,'color',color,'linewidth',3);

% ice2ra
color = [1,0,0];
m = IMAU_ICE.ice2ra.time >= 100;
xdata_lo =        IMAU_ICE.ice2ra.time(    m);
xdata_hi = flipud(IMAU_ICE.ice2ra.time(    m));
ydata_lo =        IMAU_ICE.ice2ra.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice2ra.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ice2ra.time( m),'ydata',IMAU_ICE.ice2ra.xGL_av( m)/1e3,'color',color,'linewidth',3);

% ice2rr
color = [1.0,0.63,0.0];
m = IMAU_ICE.ice2rr.time >= 100;
xdata_lo =        IMAU_ICE.ice2rr.time(    m);
xdata_hi = flipud(IMAU_ICE.ice2rr.time(    m));
ydata_lo =        IMAU_ICE.ice2rr.xGL_min( m);
ydata_hi = flipud(IMAU_ICE.ice2rr.xGL_max( m));
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H1.Ax(1,2),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H1.Ax(1,2),'xdata',IMAU_ICE.ice2rr.time( m),'ydata',IMAU_ICE.ice2rr.xGL_av( m)/1e3,'color',color,'linewidth',3);

% Legends
legend(H1.Ax(1,1),'ice0','ice1r','ice1ra','ice1rr','location','southwest');
legend(H1.Ax(1,2),'IMAU-ICE mean','IMAU-ICE range','C20 mean','C20 range','location','southwest');

%% Compare plumes for resolutions, sliding laws, stress balances, and subgrid schemes

wa = [400,60,400,60];
ha = [300,300];

margins_hor = [100,15,15,15,25];
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

ylim = [325,475];
ytick = 325:25:475;

set(H2.Ax(1,1),'ylim',ylim,'ytick',ytick);
set(H2.Ax(1,3),'ylim',ylim,'ytick',ytick,'yticklabels','');
set(H2.Ax(2,1),'ylim',ylim,'ytick',ytick);
set(H2.Ax(2,3),'ylim',ylim,'ytick',ytick,'yticklabels','');

set(H2.Ax(1,1),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H2.Ax(1,3),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H2.Ax(2,1),'xlim',[0,200],'xtick',0:50:200);
set(H2.Ax(2,3),'xlim',[0,200],'xtick',0:50:200);

ylabel(H2.Ax(1,1),'x_{GL} (km)')
ylabel(H2.Ax(2,1),'x_{GL} (km)')
xlabel(H2.Ax(2,1),'Time (yr)')
xlabel(H2.Ax(2,3),'Time (yr)')

title(H2.Ax(1,1),'Resolution')
title(H2.Ax(1,3),'Sliding law')
title(H2.Ax(2,1),'Stress balance')
title(H2.Ax(2,3),'Sub-grid scheme')
% title(H2.Ax(1,1),'A')
% title(H2.Ax(1,3),'B')
% title(H2.Ax(2,1),'C')
% title(H2.Ax(2,3),'D')

% Small axeses for bars
set(H2.Ax(1,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H2.Ax(    1,2).XAxis.Visible = 'off';
H2.Ax(    1,2).YAxis.Visible = 'off';
set(H2.Ax(1,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H2.Ax(    1,4).XAxis.Visible = 'off';
H2.Ax(    1,4).YAxis.Visible = 'off';
set(H2.Ax(2,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H2.Ax(    2,2).XAxis.Visible = 'off';
H2.Ax(    2,2).YAxis.Visible = 'off';
set(H2.Ax(2,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H2.Ax(    2,4).XAxis.Visible = 'off';
H2.Ax(    2,4).YAxis.Visible = 'off';

colors = lines(3);

% A: resolutions

% Empty line objects for legend
line('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'linewidth',3,'color',colors(1,:));
line('parent',H2.Ax(1,1),'xdata',[],'ydata',[],'linewidth',3,'color',colors(2,:));

for ri = 1:length(resolutions)
  rres = ['r' resolutions{ri}];
  color = colors(ri,:);
  
  % Plume
  xdata_lo =        IMAU_ICE.(rres).time;
  xdata_hi = flipud(IMAU_ICE.(rres).time);
  ydata_lo =        IMAU_ICE.(rres).xGL_min;
  ydata_hi = flipud(IMAU_ICE.(rres).xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H2.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H2.Ax(1,1),'xdata',IMAU_ICE.(rres).time,'ydata',IMAU_ICE.(rres).xGL_av/1e3,'color',color,'linewidth',3);
  
  % Bar
  w = 1/6;
  x = (ri-1)*1/3 + 1/12;
  y_lo = IMAU_ICE.(rres).xGL_min( end);
  y_av = IMAU_ICE.(rres).xGL_av(  end);
  y_hi = IMAU_ICE.(rres).xGL_max( end);
  patch('parent',H2.Ax(1,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H2.Ax(1,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
end

% Legend
legend(H2.Ax(1,1),'5 km','2 km','location','southwest')

% B: sliding laws

% Empty line objects for legend
line('parent',H2.Ax(1,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(1,:));
line('parent',H2.Ax(1,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(2,:));
line('parent',H2.Ax(1,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(3,:));

for sli = 1:length(sliding_laws)
  slid = sliding_laws{sli};
  color = colors(sli,:);

  % Plume
  xdata_lo =        IMAU_ICE.(slid).time;
  xdata_hi = flipud(IMAU_ICE.(slid).time);
  ydata_lo =        IMAU_ICE.(slid).xGL_min;
  ydata_hi = flipud(IMAU_ICE.(slid).xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H2.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H2.Ax(1,3),'xdata',IMAU_ICE.(slid).time,'ydata',IMAU_ICE.(slid).xGL_av/1e3,'color',color,'linewidth',3);
  
  % Bar
  w = 1/6;
  x = (sli-1)*1/3 + 1/12;
  y_lo = IMAU_ICE.(slid).xGL_min( end);
  y_av = IMAU_ICE.(slid).xGL_av(  end);
  y_hi = IMAU_ICE.(slid).xGL_max( end);
  patch('parent',H2.Ax(1,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H2.Ax(1,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
end

% Legend
legend(H2.Ax(1,3),'slid1','slid2','slid3','location','southwest')

% C: stress balances

% Empty line objects for legend
line('parent',H2.Ax(2,1),'xdata',[],'ydata',[],'linewidth',3,'color',colors(1,:));
line('parent',H2.Ax(2,1),'xdata',[],'ydata',[],'linewidth',3,'color',colors(2,:));

for sti = 1:length(stress_balances)
  strb = stress_balances{sti};
  color = colors(sti,:);
  
  % Plume
  xdata_lo =        IMAU_ICE.(strb).time;
  xdata_hi = flipud(IMAU_ICE.(strb).time);
  ydata_lo =        IMAU_ICE.(strb).xGL_min;
  ydata_hi = flipud(IMAU_ICE.(strb).xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H2.Ax(2,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H2.Ax(2,1),'xdata',IMAU_ICE.(strb).time,'ydata',IMAU_ICE.(strb).xGL_av/1e3,'color',color,'linewidth',3);
  
  % Bar
  w = 1/6;
  x = (sti-1)*1/3 + 1/12;
  y_lo = IMAU_ICE.(strb).xGL_min( end);
  y_av = IMAU_ICE.(strb).xGL_av(  end);
  y_hi = IMAU_ICE.(strb).xGL_max( end);
  patch('parent',H2.Ax(2,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H2.Ax(2,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
end

% Legend
legend(H2.Ax(2,1),'DIVA','SIA/SSA','location','southwest')

% D: subgrid schemes

% Empty line objects for legend
line('parent',H2.Ax(2,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(1,:));
line('parent',H2.Ax(2,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(2,:));
line('parent',H2.Ax(2,3),'xdata',[],'ydata',[],'linewidth',3,'color',colors(3,:));

for sbgi = 1:length(subgrid_schemes)
  sbg = subgrid_schemes{sbgi};
  color = colors(sbgi,:);
  
  % Plume
  xdata_lo =        IMAU_ICE.(sbg).time;
  xdata_hi = flipud(IMAU_ICE.(sbg).time);
  ydata_lo =        IMAU_ICE.(sbg).xGL_min;
  ydata_hi = flipud(IMAU_ICE.(sbg).xGL_max);
  xdata = [xdata_lo;xdata_hi];
  ydata = [ydata_lo;ydata_hi];
  patch('parent',H2.Ax(2,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
  line('parent',H2.Ax(2,3),'xdata',IMAU_ICE.(sbg).time,'ydata',IMAU_ICE.(sbg).xGL_av/1e3,'color',color,'linewidth',3);
  
  % Bar
  w = 1/6;
  x = (sbgi-1)*1/3 + 1/12;
  y_lo = IMAU_ICE.(sbg).xGL_min( end);
  y_av = IMAU_ICE.(sbg).xGL_av(  end);
  y_hi = IMAU_ICE.(sbg).xGL_max( end);
  patch('parent',H2.Ax(2,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
  line('parent',H2.Ax(2,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);
end

% Legend
legend(H2.Ax(2,3),'NMP','FCMP','PMP','location','southwest')

%% Compare plumes for "benchmark" run

wa = [400,60,400,60];
ha = [300,300];

margins_hor = [100,15,15,15,25];
margins_ver = [80,40,40];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + sum(wa);
hf = sum(margins_ver) + sum(ha);

H3.Fig = figure('position',[300,300,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    x = sum(margins_hor(1:i)) + sum(wa(1:i-1));
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*sum(ha(1:jp-1));
    H3.Ax(i,j) = axes('parent',H3.Fig,'units','pixels','position',[x,y,wa(i),ha(j)],'fontsize',24,...
      'xgrid','on','ygrid','on','xlim',[0,200]);
  end
end

H3.Ax = H3.Ax';

ylim = [325,475];
ytick = 325:25:475;

set(H3.Ax(1,1),'ylim',ylim,'ytick',ytick);
set(H3.Ax(1,3),'ylim',ylim,'ytick',ytick,'yticklabels','');
set(H3.Ax(2,1),'ylim',ylim,'ytick',ytick);
set(H3.Ax(2,3),'ylim',ylim,'ytick',ytick,'yticklabels','');

set(H3.Ax(1,1),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H3.Ax(1,3),'xlim',[0,200],'xtick',0:50:200,'xticklabels','');
set(H3.Ax(2,1),'xlim',[0,200],'xtick',0:50:200);
set(H3.Ax(2,3),'xlim',[0,200],'xtick',0:50:200);

ylabel(H3.Ax(1,1),'x_{GL} (km)')
ylabel(H3.Ax(2,1),'x_{GL} (km)')
xlabel(H3.Ax(2,1),'Time (yr)')
xlabel(H3.Ax(2,3),'Time (yr)')

title(H3.Ax(1,1),'Resolution')
title(H3.Ax(1,3),'Sliding law')
title(H3.Ax(2,1),'Stress balance')
title(H3.Ax(2,3),'Sub-grid scheme')
% title(H3.Ax(1,1),'A')
% title(H3.Ax(1,3),'B')
% title(H3.Ax(2,1),'C')
% title(H3.Ax(2,3),'D')

% Small axeses for bars
set(H3.Ax(1,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H3.Ax(    1,2).XAxis.Visible = 'off';
H3.Ax(    1,2).YAxis.Visible = 'off';
set(H3.Ax(1,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H3.Ax(    1,4).XAxis.Visible = 'off';
H3.Ax(    1,4).YAxis.Visible = 'off';
set(H3.Ax(2,2),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H3.Ax(    2,2).XAxis.Visible = 'off';
H3.Ax(    2,2).YAxis.Visible = 'off';
set(H3.Ax(2,4),'xlim',[-0.01,1.01],'ylim',ylim,'ytick',ytick,'xticklabels','','yticklabels','','ygrid','on','xgrid','off');
H3.Ax(    2,4).XAxis.Visible = 'off';
H3.Ax(    2,4).YAxis.Visible = 'off';

% Empty objects for legend
patch('parent',H3.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(1,:),'facealpha',0.45,'edgecolor','none');
patch('parent',H3.Ax(1,1),'xdata',[],'ydata',[],'facecolor',colors(2,:),'facealpha',0.45,'edgecolor','none');

% 5 km baseline
color = colors(1,:);

% Resolution

% Plume
xdata_lo =        benchmark_5km.rres.time;
xdata_hi = flipud(benchmark_5km.rres.time);
ydata_lo =        benchmark_5km.rres.xGL_min;
ydata_hi = flipud(benchmark_5km.rres.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(1,1),'xdata',benchmark_5km.rres.time,'ydata',benchmark_5km.rres.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 0/2 + 1/8;
y_lo = benchmark_5km.rres.xGL_min( end);
y_av = benchmark_5km.rres.xGL(     end);
y_hi = benchmark_5km.rres.xGL_max( end);
patch('parent',H3.Ax(1,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(1,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Sliding law

% Plume
xdata_lo =        benchmark_5km.slid.time;
xdata_hi = flipud(benchmark_5km.slid.time);
ydata_lo =        benchmark_5km.slid.xGL_min;
ydata_hi = flipud(benchmark_5km.slid.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(1,3),'xdata',benchmark_5km.slid.time,'ydata',benchmark_5km.slid.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 0/2 + 1/8;
y_lo = benchmark_5km.slid.xGL_min( end);
y_av = benchmark_5km.slid.xGL(     end);
y_hi = benchmark_5km.slid.xGL_max( end);
patch('parent',H3.Ax(1,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(1,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Stress balance

% Plume
xdata_lo =        benchmark_5km.strb.time;
xdata_hi = flipud(benchmark_5km.strb.time);
ydata_lo =        benchmark_5km.strb.xGL_min;
ydata_hi = flipud(benchmark_5km.strb.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(2,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(2,1),'xdata',benchmark_5km.strb.time,'ydata',benchmark_5km.strb.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 0/2 + 1/8;
y_lo = benchmark_5km.strb.xGL_min( end);
y_av = benchmark_5km.strb.xGL(     end);
y_hi = benchmark_5km.strb.xGL_max( end);
patch('parent',H3.Ax(2,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(2,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Sub-grid scheme

% Plume
xdata_lo =        benchmark_5km.sbg.time;
xdata_hi = flipud(benchmark_5km.sbg.time);
ydata_lo =        benchmark_5km.sbg.xGL_min;
ydata_hi = flipud(benchmark_5km.sbg.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(2,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(2,3),'xdata',benchmark_5km.sbg.time,'ydata',benchmark_5km.sbg.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 0/2 + 1/8;
y_lo = benchmark_5km.sbg.xGL_min( end);
y_av = benchmark_5km.sbg.xGL(     end);
y_hi = benchmark_5km.sbg.xGL_max( end);
patch('parent',H3.Ax(2,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(2,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% 2 km baseline
color = colors(2,:);

% Resolution

% Plume
xdata_lo =        benchmark_2km.rres.time;
xdata_hi = flipud(benchmark_2km.rres.time);
ydata_lo =        benchmark_2km.rres.xGL_min;
ydata_hi = flipud(benchmark_2km.rres.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(1,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(1,1),'xdata',benchmark_2km.rres.time,'ydata',benchmark_2km.rres.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 1/2 + 1/8;
y_lo = benchmark_2km.rres.xGL_min( end);
y_av = benchmark_2km.rres.xGL(     end);
y_hi = benchmark_2km.rres.xGL_max( end);
patch('parent',H3.Ax(1,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(1,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Sliding law

% Plume
xdata_lo =        benchmark_2km.slid.time;
xdata_hi = flipud(benchmark_2km.slid.time);
ydata_lo =        benchmark_2km.slid.xGL_min;
ydata_hi = flipud(benchmark_2km.slid.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(1,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(1,3),'xdata',benchmark_2km.slid.time,'ydata',benchmark_2km.slid.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 1/2 + 1/8;
y_lo = benchmark_2km.slid.xGL_min( end);
y_av = benchmark_2km.slid.xGL(     end);
y_hi = benchmark_2km.slid.xGL_max( end);
patch('parent',H3.Ax(1,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(1,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Stress balance

% Plume
xdata_lo =        benchmark_2km.strb.time;
xdata_hi = flipud(benchmark_2km.strb.time);
ydata_lo =        benchmark_2km.strb.xGL_min;
ydata_hi = flipud(benchmark_2km.strb.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(2,1),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(2,1),'xdata',benchmark_2km.strb.time,'ydata',benchmark_2km.strb.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 1/2 + 1/8;
y_lo = benchmark_2km.strb.xGL_min( end);
y_av = benchmark_2km.strb.xGL(     end);
y_hi = benchmark_2km.strb.xGL_max( end);
patch('parent',H3.Ax(2,2),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(2,2),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Sub-grid scheme

% Plume
xdata_lo =        benchmark_2km.sbg.time;
xdata_hi = flipud(benchmark_2km.sbg.time);
ydata_lo =        benchmark_2km.sbg.xGL_min;
ydata_hi = flipud(benchmark_2km.sbg.xGL_max);
xdata = [xdata_lo;xdata_hi];
ydata = [ydata_lo;ydata_hi];
patch('parent',H3.Ax(2,3),'xdata',xdata,'ydata',ydata/1e3,'facecolor',color,'facealpha',0.25,'edgecolor','none');
line('parent',H3.Ax(2,3),'xdata',benchmark_2km.sbg.time,'ydata',benchmark_2km.sbg.xGL/1e3,'color',color,'linewidth',3);
  
% Bar
w = 1/4;
x = 1/2 + 1/8;
y_lo = benchmark_2km.sbg.xGL_min( end);
y_av = benchmark_2km.sbg.xGL(     end);
y_hi = benchmark_2km.sbg.xGL_max( end);
patch('parent',H3.Ax(2,4),'xdata',[x,x+w,x+w,x],'ydata',[y_lo,y_lo,y_hi,y_hi]/1e3,'edgecolor','none','facecolor',color,'facealpha',0.45);
line('parent',H3.Ax(2,4),'xdata',[x,x+w],'ydata',[0,0]+y_av/1e3,'linewidth',4,'color',color);

% Legend
legend(H3.Ax(1,1),'5 km','2 km','location','southwest');