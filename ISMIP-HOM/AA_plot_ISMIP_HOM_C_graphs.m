clc
clear all
close all

%% Read other model data
foldername ='/Users/berends/Documents/Datasets/ISMIP-HOM/tc-2-95-2008-supplement/ismip_all/';

model_types = {...
  'aas1','FS';...
  'aas2','FS';...
  'ahu1','HO';...
  'ahu2','HO';...
  'bds1','HO';...
  'cma1','FS';...
  'cma2','HO';...
  'dpo1','HO';...
  'fpa1','HO';...
  'fpa2','FS';...
  'fsa1','HO';...
  'ghg1','FS';...
  'jvj1','FS';...
  'lpe1','HO';...
  'mbr1','HO';...
  'mmr1','FS';...
  'mtk1','HO';...
  'oga1','FS';...
  'oso1','HO';...
  'rhi1','FS';...
  'rhi2','HO';...
  'rhi3','FS';...
  'rhi4','HO';...
  'rhi5','HO';...
  'spr1','FS';...
  'ssu1','FS';...
  'tpa1','HO';...
  'yko1','FS'};

models = dir(foldername);
models = models(3:end);

experiments.c160 = [];
experiments.c080 = [];
experiments.c040 = [];
experiments.c020 = [];
experiments.c010 = [];
experiments.c005 = [];

for mi = 1: length(models)
  
  modeldata = dir([foldername models(mi).name]);
  modeldata = modeldata(3:end);
  
  % Go over all experiments, check if this model has them.
  flds = fields(experiments);
  for xi = 1:length(flds)
    ex = flds{xi};
    
    for di = 1:length(modeldata)
      mdname = modeldata(di).name;
      str = [ex '.txt'];
      if length(mdname) >= length(str)
        if strcmpi(mdname(end-length(str)+1:end),str)
          % This is the experiment from this model
          
          disp(['Reading data from model ' models(mi).name ', experiment ' ex])
          
          fid = fopen([foldername models(mi).name '/' mdname]);
          temp = textscan(fid,'%s','delimiter','\n','MultipleDelimsAsOne',1); temp = temp{1};
          fclose(fid);
          
          n = length(temp);
          nx = sqrt(n);
          if nx-floor(nx)>0
            error('whaa!')
          end
          x_vec = zeros(n,1);
          y_vec = zeros(n,1);
          u_vec = zeros(n,1);
          v_vec = zeros(n,1);
          w_vec = zeros(n,1);
          for i = 1:n
            temp2 = textscan(temp{i},'%f %f %f %f %f %f %f %f %f %f');
            x_vec(i) = temp2{1};
            y_vec(i) = temp2{2};
            u_vec(i) = temp2{3};
            v_vec(i) = temp2{4};
            w_vec(i) = temp2{5};
          end
          
          experiments.(ex).(models(mi).name).x = reshape(x_vec,[nx,nx]);
          experiments.(ex).(models(mi).name).y = reshape(y_vec,[nx,nx]);
          experiments.(ex).(models(mi).name).u = reshape(u_vec,[nx,nx]);
          experiments.(ex).(models(mi).name).v = reshape(v_vec,[nx,nx]);
          experiments.(ex).(models(mi).name).w = reshape(w_vec,[nx,nx]);
          
          if (experiments.(ex).(models(mi).name).x(1,1) == experiments.(ex).(models(mi).name).x(end,1))
            experiments.(ex).(models(mi).name).x = experiments.(ex).(models(mi).name).x';
            experiments.(ex).(models(mi).name).y = experiments.(ex).(models(mi).name).y';
            experiments.(ex).(models(mi).name).u = experiments.(ex).(models(mi).name).u';
            experiments.(ex).(models(mi).name).v = experiments.(ex).(models(mi).name).v';
            experiments.(ex).(models(mi).name).w = experiments.(ex).(models(mi).name).w';
          end
          
        end
      end
    end
  end
  
end

%% Read my own results
experiments.c160.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_160_DIVA';
experiments.c080.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_080_DIVA';
experiments.c040.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_040_DIVA';
experiments.c020.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_020_DIVA';
experiments.c010.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_010_DIVA';
experiments.c005.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_005_DIVA';

experiments.c160.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_160_hybrid';
experiments.c080.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_080_hybrid';
experiments.c040.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_040_hybrid';
experiments.c020.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_020_hybrid';
experiments.c010.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_010_hybrid';
experiments.c005.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_C_005_hybrid';

flds = fields(experiments);
for xi = 1:length(flds)
  ex = flds{xi};
  
  filename = [experiments.(ex).IMAU_ICE_DIVA.foldername '/help_fields_ANT.nc'];
  if ~exist(filename,'file')
    experiments.(ex).IMAU_ICE_DIVA.x = [];
    experiments.(ex).IMAU_ICE_DIVA.y = [];
    experiments.(ex).IMAU_ICE_DIVA.u = [];
    experiments.(ex).IMAU_ICE_DIVA.v = [];
    experiments.(ex).IMAU_ICE_DIVA.w = [];
  else
    time = ncread(filename,'time');
    ti = length(time);

    experiments.(ex).IMAU_ICE_DIVA.x = ncread(filename,'x');
    experiments.(ex).IMAU_ICE_DIVA.y = ncread(filename,'y');
    experiments.(ex).IMAU_ICE_DIVA.u = ncread(filename,'u_3D',[1,1,1,ti],[Inf,Inf,1,1]);
    experiments.(ex).IMAU_ICE_DIVA.v = ncread(filename,'v_3D',[1,1,1,ti],[Inf,Inf,1,1]);
    experiments.(ex).IMAU_ICE_DIVA.w = ncread(filename,'w_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  end
  
  filename = [experiments.(ex).IMAU_ICE_hybrid.foldername '/help_fields_ANT.nc'];
  if ~exist(filename,'file')
    experiments.(ex).IMAU_ICE_hybrid.x = [];
    experiments.(ex).IMAU_ICE_hybrid.y = [];
    experiments.(ex).IMAU_ICE_hybrid.u = [];
    experiments.(ex).IMAU_ICE_hybrid.v = [];
    experiments.(ex).IMAU_ICE_hybrid.w = [];
  else

    time = ncread(filename,'time');
    ti = length(time);

    experiments.(ex).IMAU_ICE_hybrid.x = ncread(filename,'x');
    experiments.(ex).IMAU_ICE_hybrid.y = ncread(filename,'y');
    experiments.(ex).IMAU_ICE_hybrid.u = ncread(filename,'u_3D',[1,1,1,ti],[Inf,Inf,1,1]);
    experiments.(ex).IMAU_ICE_hybrid.v = ncread(filename,'v_3D',[1,1,1,ti],[Inf,Inf,1,1]);
    experiments.(ex).IMAU_ICE_hybrid.w = ncread(filename,'w_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  end
  
end

% save('tempdata.mat','experiments')
% load('tempdata.mat')

%% Plot results
wa = 400;
ha = 300;
margin_leftleft = 50;
margin_left = 50;
margin_right = 25;
margin_botbot = 70;
margin_bot = 10;
margin_top = 50;
wf = margin_leftleft + (margin_left + wa + margin_right) * 3;
hf = margin_botbot   + (margin_bot + ha + margin_top) * 2;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
x1 = margin_leftleft + margin_left;
x2 = margin_leftleft + margin_left + wa + margin_right + margin_left;
x3 = margin_leftleft + margin_left + wa + margin_right + margin_left + wa + margin_right + margin_left;
y1 = margin_botbot + margin_bot + ha + margin_top + margin_bot;
y2 = margin_botbot + margin_bot;

H.Ax(1) = axes('units','pixels','position',[x1,y1,wa,ha],'xticklabels',[]);
H.Ax(2) = axes('units','pixels','position',[x2,y1,wa,ha],'xticklabels',[]);
H.Ax(3) = axes('units','pixels','position',[x3,y1,wa,ha],'xticklabels',[]);
H.Ax(4) = axes('units','pixels','position',[x1,y2,wa,ha]);
H.Ax(5) = axes('units','pixels','position',[x2,y2,wa,ha]);
H.Ax(6) = axes('units','pixels','position',[x3,y2,wa,ha]);

for a = 1:6
  set(H.Ax(a),'fontsize',24,'xlim',[0,1],'xgrid','on','ygrid','on');
end

title(H.Ax(1),'160 km')
title(H.Ax(2),'80 km')
title(H.Ax(3),'40 km')
title(H.Ax(4),'20 km')
title(H.Ax(5),'10 km')
title(H.Ax(6),'5 km')

ylabel(H.Ax(1),'Velocity (m/yr)')
ylabel(H.Ax(4),'Velocity (m/yr)')

xlabel(H.Ax(4),'x / L')
xlabel(H.Ax(5),'x / L')
xlabel(H.Ax(6),'x / L')

set(H.Ax(1),'ylim',[0,175]);
set(H.Ax(2),'ylim',[0,80]);
set(H.Ax(3),'ylim',[0,40]);
set(H.Ax(4),'ylim',[0,20]);
set(H.Ax(5),'ylim',[0,20]);
set(H.Ax(6),'ylim',[0,20]);

% Legend
c_IMAU = [1.0,0.2,0.2];
c_FS   = [0.2,0.5,1.0];
c_HO   = [0.1,0.7,0.3];
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_IMAU,'linewidth',3,'linestyle','--');
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_IMAU,'linewidth',3,'linestyle','-');
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7)
patch(H.Ax(1),'vertices',[],'faces',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.7)
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
line(H.Ax(1),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);

% Results
for a = 1:6
  
  if a == 1
    ex = 'c160';
  elseif a == 2
    ex = 'c080';
  elseif a == 3
    ex = 'c040';
  elseif a == 4
    ex = 'c020';
  elseif a == 5
    ex = 'c010';
  elseif a == 6
    ex = 'c005';
  end
  
  x_FS = linspace(0,1,200)';
  u_FS = [];
  u_HO = [];
  
  patch_HO = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.5);
  patch_FS = patch(H.Ax(a),'xdata',[],'ydata',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7);
  line_HO  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);
  line_FS  = line('parent',H.Ax(a),'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
  
  flds = fields(experiments.(ex));
  for mi = 1:length(flds)
    m = flds{mi};
    
    % Get transect at y = 0.25
    if strcmpi(m,'IMAU_ICE_DIVA') || strcmpi(m,'IMAU_ICE_hybrid')
      x = experiments.(ex).(m).x;
      if isempty(x); continue; end
      y = experiments.(ex).(m).y;
      yi = round(0.25*length(y));
      u = experiments.(ex).(m).u(:,yi);
      x = (x - min(x)) / range(x);
      if strcmpi(m,'IMAU_ICE_DIVA')
        ls = '-';
      else
        ls = '--';
%         % ======
%         % == hybrid solver can't handle this experiment well, skip it
%         % ======
%         continue
      end
      line('parent',H.Ax(a),'xdata',x,'ydata',u,'color',c_IMAU,'linestyle',ls,'linewidth',3);
      
    else
      x = experiments.(ex).(m).x(:,1);
      y = experiments.(ex).(m).y(1,:)';
      yi = round(0.25*length(y));
      u = experiments.(ex).(m).u(:,yi);
      
      % Determine if this model is FS or HO
      FS = false;
      HO = false;
      for mii = 1:size(model_types,1)
        if strcmpi(model_types{mii,1},m)
          if strcmpi(model_types{mii,2},'FS')
            FS = true;
          else
            HO = true;
          end
        end
      end
      if ~(FS || HO)
        for mii = 1:size(model_types,1)
          if strcmpi(model_types{mii,1}(1:3),m(1:3))
            if strcmpi(model_types{mii,2},'FS')
              FS = true;
            else
              HO = true;
            end
          end
        end
      end
      if ~(FS || HO)
        % Unknown model?
        continue
      end
      
      % Add to data ranges for HO/FS models
      up = interp1(x,u,x_FS);
      if FS
        u_FS(:,end+1) = up;
      else
        u_HO(:,end+1) = up;
      end
      
    end
    
  end
  
  % Missing data points
  m = true(size(x_FS));
  for i = 1:length(x_FS)
    if sum(isnan(u_FS(i,:)))+sum(isnan(u_HO(i,:)))>0
      m(i) = false;
    end
  end
  u_FS = u_FS(m,:);
  u_HO = u_HO(m,:);
  x_FS = x_FS(m);
  
  % ISMIP-HOM ensemble data
  uav_FS = mean(u_FS,2);
  uav_HO = mean(u_HO,2);
  sigma_FS = zeros(size(x_FS));
  sigma_HO = zeros(size(x_FS));
  for i = 1:size(u_FS,1)
    sigma_FS(i) = std(u_FS(i,:));
    sigma_HO(i) = std(u_HO(i,:));
    if isnan(sigma_FS(i)); sigma_FS(i) = 0; end
    if isnan(sigma_HO(i)); sigma_HO(i) = 0; end
  end
  umin_FS = uav_FS - sigma_FS;
  umax_FS = uav_FS + sigma_FS;
  umin_HO = uav_HO - sigma_HO;
  umax_HO = uav_HO + sigma_HO;
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_FS;flipud(umax_FS)];
  set(patch_FS,'xdata',xdata,'ydata',ydata)
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_HO;flipud(umax_HO)];
  set(patch_HO,'xdata',xdata,'ydata',ydata)
  
  set(line_HO,'xdata',x_FS,'ydata',uav_HO)
  set(line_FS,'xdata',x_FS,'ydata',uav_FS)
  
end

legend(H.Ax(1),'SIA/SSA','DIVA','Full-Stokes','Higher-Order','FS mean','HO mean','location','northwest')
% legend(H.Ax(1),'DIVA','Full-Stokes','Higher-Order','FS mean','HO mean','location','northwest')