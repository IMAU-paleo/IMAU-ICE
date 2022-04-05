clc
clear all
close all

foldernames = {...
  'exp_II_retreat_2km_perfect',...
  'exp_II_retreat_2km_visc_hi',...
  'exp_II_retreat_2km_visc_lo',...
  'exp_II_retreat_2km_SMB_hi',...
  'exp_II_retreat_2km_SMB_lo'};

%% Read data
for fi = 1:length(foldernames)
  
  filename = [foldernames{fi} '/MISMIPplus_output.nc'];

  time = ncread(filename,'time');
  xGL  = ncread(filename,'xGL');
  ny  = size(xGL,2);
  jmid = ceil(ny/2);
  xGL = xGL(:,jmid);

  results(fi).time = time;
  results(fi).xGL  = xGL;
end

%% Set up GUI

wa = 600;
ha = 400;

margins_hor = [100,25];
margins_ver = [80,25];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[100,100,wf,hf],'color','w');

for i = 1:nax
  for j = 1:nay
    x = sum(margins_hor(1:i)) + (i-1)*wa;
    jp = nay+1-j;
    y = sum(margins_ver(1:jp)) + (jp-1)*ha;
    
    ax = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24,...
      'xgrid','on','ygrid','on');
    
    H.Ax(j,i) = ax;
    
  end
end

xlabel(H.Ax,'Time (yr)')
ylabel(H.Ax,'x_{GL} (km)')

%% Plot results
colors = lines(2);
line('parent',H.Ax,'xdata',results(1).time,'ydata',results(1).xGL/1e3,'color','k'        ,'linewidth',3,'linestyle','-');
line('parent',H.Ax,'xdata',results(2).time,'ydata',results(2).xGL/1e3,'color',colors(1,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results(3).time,'ydata',results(3).xGL/1e3,'color',colors(1,:),'linewidth',3,'linestyle',':');
line('parent',H.Ax,'xdata',results(4).time,'ydata',results(4).xGL/1e3,'color',colors(2,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results(5).time,'ydata',results(5).xGL/1e3,'color',colors(2,:),'linewidth',3,'linestyle',':');

legend('Target','Visc hi','Visc lo','SMB hi','SMB lo','location','northeast')