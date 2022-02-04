clc
clear all
close all

foldernames = {...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_5km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_4km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_3km'};

% Schoof 2006 solution
n_flow = 3;
sec_per_year = 31556943.36;
A = (3.7E8 ^ (-n_flow)) * sec_per_year;
B = A^(-1/n_flow);
tantheta = 0.001;
H0 = 2000;
L = 40000;
m = 1;
y = ncread([foldernames{end} '/help_fields_ANT.nc'],'y');
[u_analytical,~] = Schoof2006_SSA_solution( 0.001, 2000, B, L, m, y);

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

colors = lines(length(foldernames));
% Empty objects for legend
line('parent',H.Ax,'xdata',[],'ydata',[],'color','k','linewidth',2);
for fi = 1:length(foldernames)
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  line('parent',H.Ax,'xdata',[],'ydata',[],'color',colors(fi,:),'linewidth',2);
end

line('parent',H.Ax,'xdata',y/1e3,'ydata',u_analytical,'linewidth',2,'color','k');
set(H.Ax,'xlim',[min(y),max(y)]/1e3);
xlabel(H.Ax,'Distance from centre (km)');
ylabel(H.Ax,'Ice velocity (m/yr)');

for fi = 1:length(foldernames)
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  
  time = ncread(filename,'time'); ti = length(time);
  x    = ncread(filename,'x');    xi = round(length(x)/2);
  y    = ncread(filename,'y');
  u    = ncread(filename,'u_base',[xi,1,ti],[1,Inf,1]);
  
  line('parent',H.Ax,'xdata',y/1e3,'ydata',u,'color',colors(fi,:),'linewidth',2);
end

legend(H.Ax,'Schoof2006','dx = 5 km','dx = 4 km','dx = 3km');
% legend(H.Ax,'Schoof2006','dx = 5 km','dx = 4 km','dx = 3km','dx = 2km','dx = 1km','dx = 500 m','dx = 250 m');