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

experiments.e000 = [];
experiments.e001 = [];

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
          x_vec = zeros(n,1);
          u_vec = zeros(n,1);
          w_vec = zeros(n,1);
          
          for i = 1:n
            temp2 = textscan(temp{i},'%f %f %f %f %f');
            x_vec(i) = temp2{1};
            u_vec(i) = temp2{2};
            w_vec(i) = temp2{3};
          end
          
          u_vec = flipud(u_vec);
          w_vec = flipud(w_vec);
          
          if max(x_vec)>100
            x_vec = x_vec/max(x_vec);
          end
          
          experiments.(ex).(models(mi).name).x = x_vec;
          experiments.(ex).(models(mi).name).u = u_vec;
          experiments.(ex).(models(mi).name).w = w_vec;
          
%           line(x_vec,u_vec);
          
        end
      end
    end
  end
  
end

% Read my own results
experiments.e000.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_E_0_DIVA';
experiments.e001.IMAU_ICE_DIVA.foldername   = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_E_1_DIVA';

experiments.e000.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_E_0_hybrid';
experiments.e001.IMAU_ICE_hybrid.foldername = '/Users/berends/Documents/Models/IMAU-ICE/ISMIP-HOM/ISMIP_HOM_E_1_hybrid';

flds = fields(experiments);
for xi = 1:length(flds)
  ex = flds{xi};
  
  filename = [experiments.(ex).IMAU_ICE_DIVA.foldername '/help_fields_ANT.nc'];
  
  time = ncread(filename,'time');
  ti = length(time);
  
  experiments.(ex).IMAU_ICE_DIVA.x = ncread(filename,'x');
  experiments.(ex).IMAU_ICE_DIVA.y = ncread(filename,'y');
  experiments.(ex).IMAU_ICE_DIVA.u = ncread(filename,'u_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  experiments.(ex).IMAU_ICE_DIVA.v = ncread(filename,'v_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  experiments.(ex).IMAU_ICE_DIVA.w = ncread(filename,'w_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  
  filename = [experiments.(ex).IMAU_ICE_hybrid.foldername '/help_fields_ANT.nc'];
  
  time = ncread(filename,'time');
  ti = length(time);
  
  experiments.(ex).IMAU_ICE_hybrid.x = ncread(filename,'x');
  experiments.(ex).IMAU_ICE_hybrid.y = ncread(filename,'y');
  experiments.(ex).IMAU_ICE_hybrid.u = ncread(filename,'u_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  experiments.(ex).IMAU_ICE_hybrid.v = ncread(filename,'v_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  experiments.(ex).IMAU_ICE_hybrid.w = ncread(filename,'w_3D',[1,1,1,ti],[Inf,Inf,1,1]);
  
end

%% Plot results
wa = 400;
ha = 300;
margin_left  = 90;
margin_right = 90;
margin_mid   = 25;
margin_bot = 80;
margin_top = 50;
wf = margin_left + wa + margin_mid + wa + margin_right;
hf = margin_bot  + ha + margin_top;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
x1 = margin_left;
x2 = margin_left + wa + margin_mid;
y1 = margin_bot;

H.Ax(1) = axes('units','pixels','position',[x1,y1,wa,ha]);
H.Ax(2) = axes('units','pixels','position',[x2,y1,wa,ha]);

for a = 1:2
  set(H.Ax(a),'fontsize',24,'xlim',[0,1],'xgrid','on','ygrid','on');
end

title(H.Ax(1),'No-slip bed')
title(H.Ax(2),'Slip bed')

ylabel(H.Ax(1),'Velocity (m/yr)')
ylabel(H.Ax(2),'Velocity (m/yr)'); set(H.Ax(2),'yaxislocation','right')

xlabel(H.Ax(1),'x / L')
xlabel(H.Ax(2),'x / L')

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
for a = 1:2
  
  if a == 1
    ex = 'e000';
  elseif a == 2
    ex = 'e001';
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
    
    % Get transect at y = 0.5
    if strcmpi(m,'IMAU_ICE_DIVA') || strcmpi(m,'IMAU_ICE_hybrid')
      x = experiments.(ex).(m).x;
      y = experiments.(ex).(m).y;
      yi = round(0.5*length(y));
      u = flipud(experiments.(ex).(m).u(:,yi));
      x = (x - min(x)) / range(x);
      if strcmpi(m,'IMAU_ICE_DIVA')
        ls = '-';
      else
        ls = '--';
        ls = '--';
%         % ======
%         % == hybrid solver can't handle this experiment well, skip it
%         % ======
%         continue
      end
      line('parent',H.Ax(a),'xdata',x,'ydata',u,'color',c_IMAU,'linestyle',ls,'linewidth',3);
    else
      x = experiments.(ex).(m).x;
      u = experiments.(ex).(m).u;
      
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