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

%% Set up GUI

wa          = 400;
margins_hor = [100, 40, 25];
ha          = 300;
margins_ver = [80, 40, 40];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[100,100,wf,hf],'color','w');
for i = 1:nax
  for j = 1:nay
    jp = nay+1-j;
    x = sum(margins_hor(1:i )) + (i -1)*wa;
    y = sum(margins_ver(1:jp)) + (jp-1)*ha;
    H.Ax(j,i) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24);
  end
end

xlim  = [0,200];
xtick = 0:50:200;
ylim  = [300,475];
ytick = 300:25:475;

set(H.Ax(1,1),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(1,2),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(2,1),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');
set(H.Ax(2,2),'xlim',xlim,'ylim',ylim,'xtick',xtick,'ytick',ytick,'xgrid','on','ygrid','on');

set(H.Ax(1,1),'xticklabels',{});
set(H.Ax(1,2),'xticklabels',{});
set(H.Ax(1,2),'yticklabels',{});
set(H.Ax(2,2),'yticklabels',{});

title(H.Ax(1,1),'A');
title(H.Ax(1,2),'B');
title(H.Ax(2,1),'C');
title(H.Ax(2,2),'D');

xlabel(H.Ax(2,1),'Time (yr)')
xlabel(H.Ax(2,2),'Time (yr)')
ylabel(H.Ax(1,1),'x_{GL} (km)')
ylabel(H.Ax(2,1),'x_{GL} (km)')

%% Plot results

colors = lines(3);

% Resolution
ax = H.Ax(1,1);
R1 = [];
R2 = [];
for i = 1:IMAU_ICE.nruns
  if     strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_5km_Schoof2005_DIVA_FCMP')
    R1 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_FCMP')
    R2 = IMAU_ICE.results{i};
  end
end
line('parent',ax,'xdata',R1.time,'ydata',R1.xGL/1e3,'color',colors(1,:),'linewidth',3);
line('parent',ax,'xdata',R2.time,'ydata',R2.xGL/1e3,'color',colors(2,:),'linewidth',3);
legend(ax,'5 km','2 km','location','northeast')

% Sliding law
ax = H.Ax(1,2);
R1 = [];
R2 = [];
R3 = [];
for i = 1:IMAU_ICE.nruns
  if     strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Weertman_DIVA_FCMP')
    R1 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Tsai2015_DIVA_FCMP')
    R2 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_FCMP')
    R3 = IMAU_ICE.results{i};
  end
end
line('parent',ax,'xdata',R1.time,'ydata',R1.xGL/1e3,'color',colors(1,:),'linewidth',3);
line('parent',ax,'xdata',R2.time,'ydata',R2.xGL/1e3,'color',colors(2,:),'linewidth',3);
line('parent',ax,'xdata',R3.time,'ydata',R3.xGL/1e3,'color',colors(3,:),'linewidth',3);
legend(ax,'Weertman','Tsai2015','Schoof2005','location','northeast')

% Stress balance
ax = H.Ax(2,1);
R1 = [];
R2 = [];
for i = 1:IMAU_ICE.nruns
  if     strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_5km_Schoof2005_hybrid_FCMP')
    R1 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_FCMP')
    R2 = IMAU_ICE.results{i};
  end
end
line('parent',ax,'xdata',R1.time,'ydata',R1.xGL/1e3,'color',colors(1,:),'linewidth',3);
line('parent',ax,'xdata',R2.time,'ydata',R2.xGL/1e3,'color',colors(2,:),'linewidth',3);
legend(ax,'DIVA','hybrid','location','northeast')

% Subgrid scheme
ax = H.Ax(2,2);
% Empty line objects for legend
line('parent',ax,'xdata',[],'ydata',[],'color',colors(1,:),'linewidth',3)
line('parent',ax,'xdata',[],'ydata',[],'color',colors(2,:),'linewidth',3)
line('parent',ax,'xdata',[],'ydata',[],'color',colors(3,:),'linewidth',3)
line('parent',ax,'xdata',[],'ydata',[],'color','k','linewidth',3,'linestyle','--')
line('parent',ax,'xdata',[],'ydata',[],'color','k','linewidth',3,'linestyle','-')
R1 = [];
R2 = [];
R3 = [];
R4 = [];
R5 = [];
R6 = [];
for i = 1:IMAU_ICE.nruns
  if     strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_5km_Schoof2005_DIVA_NMP')
    R1 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_5km_Schoof2005_DIVA_FCMP')
    R2 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_5km_Schoof2005_DIVA_PMP')
    R3 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_NMP')
    R4 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_FCMP')
    R5 = IMAU_ICE.results{i};
  elseif strcmpi(IMAU_ICE.results{i}.run,'MISMIPplus_ice1r_2km_Schoof2005_DIVA_PMP')
    R6 = IMAU_ICE.results{i};
  end
end
line('parent',ax,'xdata',R1.time,'ydata',R1.xGL/1e3,'color',colors(1,:),'linewidth',3,'linestyle','--');
line('parent',ax,'xdata',R2.time,'ydata',R2.xGL/1e3,'color',colors(2,:),'linewidth',3,'linestyle','--');
line('parent',ax,'xdata',R3.time,'ydata',R3.xGL/1e3,'color',colors(3,:),'linewidth',3,'linestyle','--');
line('parent',ax,'xdata',R4.time,'ydata',R4.xGL/1e3,'color',colors(1,:),'linewidth',3);
line('parent',ax,'xdata',R5.time,'ydata',R5.xGL/1e3,'color',colors(2,:),'linewidth',3);
line('parent',ax,'xdata',R6.time,'ydata',R6.xGL/1e3,'color',colors(3,:),'linewidth',3);
legend(ax,'NMP','FCMP','PMP','5 km','2 km','location','northeast')
