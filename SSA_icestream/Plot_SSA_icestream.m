clc
clear all
close all

foldernames = {...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_40km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_32km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_20km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_16km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_10km'};

% Schoof 2006 solution
n_flow = 3;
sec_per_year = 31556943.36;
A = 1e-18;
B = A^(-1/n_flow);
tantheta = 0.0003;
H0 = 2000;
L = 150000;
m = 1;
y = ncread([foldernames{end} '/help_fields_ANT.nc'],'y');
[u_analytical,~] = Schoof2006_SSA_solution( tantheta, H0, B, L, m, y);

%% Plot

ha = 500;
wa = 700;
margin_left  = 100;
margin_right = 25;
margin_bot   = 100;
margin_top   = 50;

wf = margin_left + wa + margin_right;
hf = margin_bot + ha + margin_top;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('units','pixels','position',[margin_left,margin_bot,wa,ha],...
  'fontsize',24,'xgrid','on','ygrid','on');

xlabel(H.Ax,'Distance from centre (km)');
ylabel(H.Ax,'Ice velocity (m/yr)');
legend(H.Ax,'Schoof (2006) analytical solution')

% Empty objects for legend
colors = lines(length(foldernames));

line('parent',H.Ax,'xdata',[],'ydata',[],'color','k','linewidth',2);
for fi = 1:length(foldernames)
  line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(fi,:),'linewidth',2);
end

% Analytical solution
line('parent',H.Ax,'xdata',y/1e3,'ydata',u_analytical,'color','k','linewidth',2);
set(H.Ax,'xlim',[min(y),max(y)]/1e3);

% Modelled solutions
for fi = 1:length(foldernames)
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  
  time = ncread(filename,'time'); ti = length(time);
  x    = ncread(filename,'x');    xi = round(length(x)/2);
  y    = ncread(filename,'y');
  u    = ncread(filename,'u_surf',[xi,1,ti],[1,Inf,1]);
  
  line('parent',H.Ax,'xdata',y/1e3,'ydata',u,'color',colors(fi,:),'linewidth',2);
end

legend(H.Ax,'Schoof 2006','dx = 40 km','dx = 32 km','dx = 20 km','dx = 16 km','dx = 10 km');