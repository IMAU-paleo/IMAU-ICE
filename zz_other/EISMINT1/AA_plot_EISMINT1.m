clc
clear all
close all

resstring = '50km';

ensemble_name = 'EISMINT1';
experiments = {'A','B','C','D','E','F'};

%% Read IMAU-ICE results
for xpi = 1:length(experiments)
  xp = experiments{xpi};
  foldername = [ensemble_name '_' xp '_' resstring];
  filename_restart     = [foldername '/restart_ANT.nc'];
  filename_help_fields = [foldername '/help_fields_ANT.nc'];
  
  % Read dimensions
  time = ncread(filename_restart,'time');
  x    = ncread(filename_restart,'x');
  y    = ncread(filename_restart,'y');
  zeta = ncread(filename_restart,'zeta');
  
  % Get indices of ice divide grid cell
  i_mid = find(x==0);
  j_mid = find(y==0);
  nz = length(zeta);
  
  % Read model data
  results.(xp).time        = time;
  results.(xp).Hi          = permute(ncread(filename_restart,    'Hi'    ,[i_mid,j_mid,   1],[1,1,  Inf]),[3,1,2  ]);
  results.(xp).Ti_base     = permute(ncread(filename_help_fields,'Ti'    ,[i_mid,j_mid,nz,1],[1,1,1,Inf]),[4,1,2,3]);
  results.(xp).Ti_pmp_base = permute(ncread(filename_help_fields,'Ti_pmp',[i_mid,j_mid,nz,1],[1,1,1,Inf]),[4,1,2,3]);
  results.(xp).Ti_base_rel = results.(xp).Ti_base - results.(xp).Ti_pmp_base;

end

%% Ice thickness change at divide
wa = 600;
ha = 500;

margins_hor = [100,50,25];
margins_ver = [100,50];

wf = sum(margins_hor) + (length(margins_hor)-1)*wa;
hf = sum(margins_ver) + (length(margins_ver)-1)*ha;

H1.Fig = figure('position',[200,200,wf,hf],'color','w');
H1.Ax(1) = axes('parent',H1.Fig,'units','pixels','position',[sum(margins_hor(1:1))+wa*0, sum(margins_ver(1:1))+ha*0, wa, ha]);
H1.Ax(2) = axes('parent',H1.Fig,'units','pixels','position',[sum(margins_hor(1:2))+wa*1, sum(margins_ver(1:1))+ha*0, wa, ha]);

for a = 1:2
  set(H1.Ax(a),'fontsize',24,'xlim',[0,120],'xtick',0:20:120,'ylim',[-400,400],'ytick',-400:100:400,'xgrid','on','ygrid','on');
end
set(H1.Ax(2),'yticklabels',{});

xlabel(H1.Ax(1),'Time (kyr)')
xlabel(H1.Ax(2),'Time (kyr)')
ylabel(H1.Ax(1),'Ice thickness change at divide (m)');
title(H1.Ax(1),'a');
title(H1.Ax(2),'b');

% empty line objects for legend
line('parent',H1.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','b');
line('parent',H1.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','r');
line('parent',H1.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','b');
line('parent',H1.Ax(2),'xdata',[],'ydata',[],'linewidth',2,'color','r');

% plot data
line('parent',H1.Ax(1),'xdata',results.B.time/1e3,'ydata',results.B.Hi - results.B.Hi(1),'linewidth',2,'color','b');
line('parent',H1.Ax(1),'xdata',results.C.time/1e3,'ydata',results.C.Hi - results.C.Hi(1),'linewidth',2,'color','r');
line('parent',H1.Ax(2),'xdata',results.E.time/1e3,'ydata',results.E.Hi - results.E.Hi(1),'linewidth',2,'color','b');
line('parent',H1.Ax(2),'xdata',results.F.time/1e3,'ydata',results.F.Hi - results.F.Hi(1),'linewidth',2,'color','r');

% find G-IG ranges
m_20kyr = results.B.time > 100000;
m_40kyr = results.B.time >  80000;
A_Hi_moving_20kyr = range( results.B.Hi( m_20kyr));
A_Hi_moving_40kyr = range( results.C.Hi( m_40kyr));
A_Hi_fixed_20kyr  = range( results.E.Hi( m_20kyr));
A_Hi_fixed_40kyr  = range( results.F.Hi( m_40kyr));

% legend with Huybrechts ranges
legend(H1.Ax(1),...
  ['20 kyr: R = ' num2str(round(A_Hi_moving_20kyr*10)/10) ' m (512.0 - 551.5 m)'],...
  ['40 kyr: R = ' num2str(round(A_Hi_moving_40kyr*10)/10) ' m (574.9 - 616.0 m)']);
legend(H1.Ax(2),...
  ['20 kyr: R = ' num2str(round(A_Hi_fixed_20kyr *10)/10) ' m (552.5 - 585.3 m)'],...
  ['40 kyr: R = ' num2str(round(A_Hi_fixed_40kyr *10)/10) ' m (614.2 - 649.5 m)']);

%% Relative basal temperature change at divide
wa = 600;
ha = 500;

margins_hor = [100,50,25];
margins_ver = [100,50];

wf = sum(margins_hor) + (length(margins_hor)-1)*wa;
hf = sum(margins_ver) + (length(margins_ver)-1)*ha;

H2.Fig = figure('position',[400,400,wf,hf],'color','w');
H2.Ax(1) = axes('parent',H2.Fig,'units','pixels','position',[sum(margins_hor(1:1))+wa*0, sum(margins_ver(1:1))+ha*0, wa, ha]);
H2.Ax(2) = axes('parent',H2.Fig,'units','pixels','position',[sum(margins_hor(1:2))+wa*1, sum(margins_ver(1:1))+ha*0, wa, ha]);

for a = 1:2
  set(H2.Ax(a),'fontsize',24,'xlim',[0,120],'xtick',0:20:120,'ylim',[-4,6],'ytick',-4:2:6,'xgrid','on','ygrid','on');
end
set(H2.Ax(2),'yticklabels',{});

xlabel(H2.Ax(1),'Time (kyr)')
xlabel(H2.Ax(2),'Time (kyr)')
ylabel(H2.Ax(1),'Relative basal temperature change at divide (K)');
title(H2.Ax(1),'a');
title(H2.Ax(2),'b');

% empty line objects for legend
line('parent',H2.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','b');
line('parent',H2.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','r');
line('parent',H2.Ax(1),'xdata',[],'ydata',[],'linewidth',2,'color','b');
line('parent',H2.Ax(2),'xdata',[],'ydata',[],'linewidth',2,'color','r');

% plot data
line('parent',H2.Ax(1),'xdata',results.B.time/1e3,'ydata',results.B.Ti_base_rel - results.B.Ti_base_rel(1),'linewidth',2,'color','b');
line('parent',H2.Ax(1),'xdata',results.C.time/1e3,'ydata',results.C.Ti_base_rel - results.C.Ti_base_rel(1),'linewidth',2,'color','r');
line('parent',H2.Ax(2),'xdata',results.E.time/1e3,'ydata',results.E.Ti_base_rel - results.E.Ti_base_rel(1),'linewidth',2,'color','b');
line('parent',H2.Ax(2),'xdata',results.F.time/1e3,'ydata',results.F.Ti_base_rel - results.F.Ti_base_rel(1),'linewidth',2,'color','r');

% find G-IG ranges
m_20kyr = results.B.time > 100000;
m_40kyr = results.B.time >  80000;
A_Ti_moving_20kyr = range( results.B.Ti_base_rel( m_20kyr));
A_Ti_moving_40kyr = range( results.C.Ti_base_rel( m_40kyr));
A_Ti_fixed_20kyr  = range( results.E.Ti_base_rel( m_20kyr));
A_Ti_fixed_40kyr  = range( results.F.Ti_base_rel( m_40kyr));

% legend with Huybrechts ranges
legend(H2.Ax(1),...
  ['20 kyr: R = ' num2str(round(A_Ti_moving_20kyr*100)/100) ' K (2.08 - 2.57 K)'],...
  ['40 kyr: R = ' num2str(round(A_Ti_moving_40kyr*100)/100) ' K (7.02 - 7.72 K)']);
legend(H2.Ax(2),...
  ['20 kyr: R = ' num2str(round(A_Ti_fixed_20kyr *100)/100) ' K (1.70 - 2.20 K)'],...
  ['40 kyr: R = ' num2str(round(A_Ti_fixed_40kyr *100)/100) ' K (2.82 - 4.18 K)']);