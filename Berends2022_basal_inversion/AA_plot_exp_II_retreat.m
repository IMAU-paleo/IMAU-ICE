clc
clear all
close all

foldernames = {...
  'exp_II_retreat_5km_target',...        %  1
  'exp_II_retreat_5km_unperturbed',...   %  2
  'exp_II_retreat_5km_visc_hi',...       %  3
  'exp_II_retreat_5km_visc_lo',...       %  4
  'exp_II_retreat_5km_SMB_hi',...        %  5
  'exp_II_retreat_5km_SMB_lo',...        %  6
  'exp_II_retreat_5km_BMB_hi',...        %  7
  'exp_II_retreat_5km_BMB_lo',...        %  8
  'exp_II_retreat_5km_topo_hi',...       %  9
  'exp_II_retreat_5km_topo_lo',...       % 10
  'exp_II_retreat_5km_ut_hi',...         % 11
  'exp_II_retreat_5km_ut_lo',...         % 12
  'exp_II_retreat_5km_p_hi',...          % 13
  'exp_II_retreat_5km_p_lo',...          % 14
  'exp_II_retreat_5km_noneq'};           % 15

%% Read data
for fi = 1:length(foldernames)
  
  % GL position
  filename = [foldernames{fi} '/MISMIPplus_output.nc'];

  time = ncread(filename,'time');
  xGL  = ncread(filename,'xGL');
  ny  = size(xGL,2);
  jmid = ceil(ny/2);
  xGL = xGL(:,jmid);

  results(fi).timeGL = time;
  results(fi).xGL    = xGL;
  
  % Ice volume
  filename = [foldernames{fi} '/restart_ANT.nc'];
  x   = ncread(filename,'x');
  dx  = x(2) - x(1);
  Hi  = ncread(filename,'Hi');
  Hb  = ncread(filename,'Hb');
  
  ice_density                      =  910.0;
  seawater_density                 = 1028.0;
  TAF = Hi - max(0, (-Hb) * (seawater_density / ice_density));
  
  results(fi).timeVi = time;
  results(fi).Vi = zeros(size(time));
  for ti = 1:length(time)
    TAFp = max(0,TAF(:,:,ti));
    results(fi).Vi(ti) = sum(sum(TAFp)) * dx^2;
  end
  
  results(fi).Vi = results(fi).Vi - results(fi).Vi(1);
%   results(fi).Vi = results(fi).Vi * 100 / results(fi).Vi(1);
  
end

%% Set up GUI

wa = 500;
ha = 350;

% margins_hor = [110,45,100];
margins_hor = [110,25];
margins_ver = [80,50];

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

% set(H.Ax(1,2),'yaxislocation','right')

% title(H.Ax(1,1),'A');
% title(H.Ax(1,2),'B');

xlabel(H.Ax(1,1),'Time (yr)')
% xlabel(H.Ax(1,2),'Time (yr)')

ylabel(H.Ax(1,1),'\Delta V_{af} (m^3)')
% ylabel(H.Ax(1,2),'x_{GL} (km)')

%% Plot results
colors = lines(7);

% Empty line objects for legend
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color','k'        ,'linewidth',4,'linestyle','-');   % Unperturbed
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color','k'        ,'linewidth',3,'linestyle','--');  % Perturbed high
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color','k'        ,'linewidth',3,'linestyle',':');   % Perturbed low
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(1,:),'linewidth',3,'linestyle','-');   % Viscosity
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(2,:),'linewidth',3,'linestyle','-');   % SMB
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(3,:),'linewidth',3,'linestyle','-');   % BMB
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(4,:),'linewidth',3,'linestyle','-');   % Topo
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(5,:),'linewidth',3,'linestyle','-');   % ut
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(6,:),'linewidth',3,'linestyle','-');   % p
line('parent',H.Ax(1,1),'xdata',[],'ydata',[],'color',colors(7,:),'linewidth',3,'linestyle','-');   % Non-eq.

% Actual results

% Viscosity
line('parent',H.Ax(1,1),'xdata',results( 3).timeVi,'ydata',results( 3).Vi     ,'color',colors(1,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results( 3).timeGL,'ydata',results( 3).xGL/1e3,'color',colors(1,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results( 4).timeVi,'ydata',results( 4).Vi     ,'color',colors(1,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results( 4).timeGL,'ydata',results( 4).xGL/1e3,'color',colors(1,:),'linewidth',3,'linestyle',':');
% SMB
line('parent',H.Ax(1,1),'xdata',results( 5).timeVi,'ydata',results( 5).Vi     ,'color',colors(2,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results( 5).timeGL,'ydata',results( 5).xGL/1e3,'color',colors(2,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results( 6).timeVi,'ydata',results( 6).Vi     ,'color',colors(2,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results( 6).timeGL,'ydata',results( 6).xGL/1e3,'color',colors(2,:),'linewidth',3,'linestyle',':');
% SMB
line('parent',H.Ax(1,1),'xdata',results( 7).timeVi,'ydata',results( 7).Vi     ,'color',colors(3,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results( 7).timeGL,'ydata',results( 7).xGL/1e3,'color',colors(3,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results( 8).timeVi,'ydata',results( 8).Vi     ,'color',colors(3,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results( 8).timeGL,'ydata',results( 8).xGL/1e3,'color',colors(3,:),'linewidth',3,'linestyle',':');
% Topo
line('parent',H.Ax(1,1),'xdata',results( 9).timeVi,'ydata',results( 9).Vi     ,'color',colors(4,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results( 9).timeGL,'ydata',results( 9).xGL/1e3,'color',colors(4,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results(10).timeVi,'ydata',results(10).Vi     ,'color',colors(4,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results(10).timeGL,'ydata',results(10).xGL/1e3,'color',colors(4,:),'linewidth',3,'linestyle',':');
% ut
line('parent',H.Ax(1,1),'xdata',results(11).timeVi,'ydata',results(11).Vi     ,'color',colors(5,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results(11).timeGL,'ydata',results(11).xGL/1e3,'color',colors(5,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results(12).timeVi,'ydata',results(12).Vi     ,'color',colors(5,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results(12).timeGL,'ydata',results(12).xGL/1e3,'color',colors(5,:),'linewidth',3,'linestyle',':');
% p
line('parent',H.Ax(1,1),'xdata',results(13).timeVi,'ydata',results(13).Vi     ,'color',colors(6,:),'linewidth',3,'linestyle','--');
% line('parent',H.Ax(1,2),'xdata',results(13).timeGL,'ydata',results(13).xGL/1e3,'color',colors(6,:),'linewidth',3,'linestyle','--');
line('parent',H.Ax(1,1),'xdata',results(14).timeVi,'ydata',results(14).Vi     ,'color',colors(6,:),'linewidth',3,'linestyle',':');
% line('parent',H.Ax(1,2),'xdata',results(14).timeGL,'ydata',results(14).xGL/1e3,'color',colors(6,:),'linewidth',3,'linestyle',':');
% Non-eq.
line('parent',H.Ax(1,1),'xdata',results(15).timeVi,'ydata',results(15).Vi     ,'color',colors(7,:),'linewidth',3,'linestyle','-');
% line('parent',H.Ax(1,2),'xdata',results(15).timeGL,'ydata',results(15).xGL/1e3,'color',colors(7,:),'linewidth',3,'linestyle','-');

% Target
line('parent',H.Ax(1,1),'xdata',results( 1).timeVi,'ydata',results( 1).Vi     ,'color','k'        ,'linewidth',4,'linestyle','-');
% line('parent',H.Ax(1,2),'xdata',results( 1).timeGL,'ydata',results( 1).xGL/1e3,'color','k'        ,'linewidth',4,'linestyle','-');

legend(H.Ax(1,1),'Target','Perturbed high','Perturbed low',...
  'Viscosity','SMB','BMB','Topo','Z-I u\_t','Z-I p','Non-eq.','location','southwest')