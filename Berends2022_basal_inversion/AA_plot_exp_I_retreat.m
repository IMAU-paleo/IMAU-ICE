clc
clear all
close all

foldernames = {...
  'exp_I_retreat_20km_target',...
  'exp_I_retreat_20km_unperturbed',...
  'exp_I_retreat_20km_visc_hi',...
  'exp_I_retreat_20km_visc_lo',...
  'exp_I_retreat_20km_SMB_hi',...
  'exp_I_retreat_20km_SMB_lo',...
  'exp_I_retreat_20km_topo_hi',...
  'exp_I_retreat_20km_topo_lo',...
  'exp_I_retreat_20km_ut_hi',...
  'exp_I_retreat_20km_ut_lo',...
  'exp_I_retreat_20km_p_hi',...
  'exp_I_retreat_20km_p_lo'};

%% Read data
ice_density                      =  910.0;
seawater_density                 = 1028.0;
ocean_area                       = 3.611E14;
  
for fi = 1:length(foldernames)
  
  filename = [foldernames{fi} '/scalar_output_ANT.nc'];

  time = ncread(filename,'time');
  Vi   = ncread(filename,'ice_volume') * seawater_density * ocean_area / ice_density;
  

  results(fi).time = time;
  results(fi).Vi   = Vi;
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
ylabel(H.Ax,'Ice volume (m^3)')

%% Plot results
colors = lines(5);

% Empty line objects for legend
line('parent',H.Ax,'xdata',[],'ydata',[],'color','k'        ,'linewidth',4,'linestyle','-');  % Unperturbed
line('parent',H.Ax,'xdata',[],'ydata',[],'color','k'        ,'linewidth',3,'linestyle','--');  % Perturbed high
line('parent',H.Ax,'xdata',[],'ydata',[],'color','k'        ,'linewidth',3,'linestyle',':');   % Perturbed low
line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(1,:),'linewidth',3,'linestyle','-');   % Viscosity
line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(2,:),'linewidth',3,'linestyle','-');   % SMB
line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(3,:),'linewidth',3,'linestyle','-');   % Topo
line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(4,:),'linewidth',3,'linestyle','-');   % ut
line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(5,:),'linewidth',3,'linestyle','-');   % p

% Actual results

% Unperturbed
line('parent',H.Ax,'xdata',results( 2).time,'ydata',results( 2).Vi,'color','k'        ,'linewidth',4,'linestyle','-');
% Viscosity
line('parent',H.Ax,'xdata',results( 3).time,'ydata',results( 3).Vi,'color',colors(1,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results( 4).time,'ydata',results( 4).Vi,'color',colors(1,:),'linewidth',3,'linestyle',':');
% SMB
line('parent',H.Ax,'xdata',results( 5).time,'ydata',results( 5).Vi,'color',colors(2,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results( 6).time,'ydata',results( 6).Vi,'color',colors(2,:),'linewidth',3,'linestyle',':');
% Topo
line('parent',H.Ax,'xdata',results( 7).time,'ydata',results( 7).Vi,'color',colors(3,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results( 8).time,'ydata',results( 8).Vi,'color',colors(3,:),'linewidth',3,'linestyle',':');
% ut
line('parent',H.Ax,'xdata',results( 9).time,'ydata',results( 9).Vi,'color',colors(4,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results(10).time,'ydata',results(10).Vi,'color',colors(4,:),'linewidth',3,'linestyle',':');
% p
line('parent',H.Ax,'xdata',results(11).time,'ydata',results(11).Vi,'color',colors(5,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax,'xdata',results(12).time,'ydata',results(12).Vi,'color',colors(5,:),'linewidth',3,'linestyle',':');

legend('Unperturbed','Perturbed high','Perturbed low','Viscosity','SMB','Topo','Z-I u\_t','Z-I p','location','northeast')