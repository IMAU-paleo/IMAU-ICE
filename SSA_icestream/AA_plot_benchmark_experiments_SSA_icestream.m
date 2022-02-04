clc
clear all
close all

foldernames = {...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_5km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_4km',...
  '/Users/berends/Documents/Models/IMAU-ICE/SSA_icestream/SSA_icestream_3km'};

%% Plot

ha  = 500;
wa1 = 500;
wa2 = 700;
margin_left  = 25;
margin_mid   = 50;
margin_right = 100;
margin_bot   = 100;
margin_top   = 50;

wf = margin_left + wa1 + margin_mid + wa2 + margin_right;
hf = margin_bot + ha + margin_top;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax1 = axes('units','pixels','position',[margin_left,margin_bot,wa1,ha]);
H.Ax1.XAxis.Visible = 'off';
H.Ax1.YAxis.Visible = 'off';
H.Ax1.ZAxis.Visible = 'off';

%% Panel 1: 3-D display of the schematic ice sheet

lx = 1;
ly = 1;
hz = 0.3;
hi = 0.6;

ar.n = 21;
ar.l = 0.2;

% ar.X = zeros(ar.n,1);
% ar.Y = zeros(ar.n,1);
% ar.Z = zeros(ar.n,1);
% ar.U = zeros(ar.n,1);
% ar.V = zeros(ar.n,1);
% ar.W = zeros(ar.n,1);

ar.yv = linspace(-120000,120000,ar.n);

j = 0;
for x0 = 0:0.25:0.75
  for i = 1:ar.n
    % Schoof 2006 solution
    n_flow = 3;
    sec_per_year = 31556943.36;
    A = (3.7E8 ^ (-n_flow)) * sec_per_year;
    B = A^(-1/n_flow);
    tantheta = 0.001;
    H0 = 2000;
    L = 40000;
    m = 1;
    y = ar.yv(i);
    [u,~] = Schoof2006_SSA_solution( 0.001, 2000, B, L, m, y);
    % Arrow coordinates
    j = j+1;
    ar.X(j) = x0;
    ar.Y(j) = (i-1)/(ar.n-1);
    ar.Z(j) = x0 * hi + (1-x0) * (hi+hz);
    ar.U(j) = (u/777.5370)*ar.l;
    ar.V(j) = 0;
    ar.W(j) = -(u/777.5370)*ar.l*(hz/lx);
  end
end

q = quiver3(H.Ax1,ar.X,ar.Y,ar.Z,ar.U,ar.V,ar.W,'color','k');

V = [0  0  0;...
     lx 0  0;...
     lx ly 0;...
     0  ly 0;...
     0  0  hz;...
     0  ly hz;...
     0  0  hz+hi;...
     1  0  hi;...
     1  1  hi;...
     0  1  hz+hi];
   
F_bed = [1 2 3 4;...
         1 2 5 1;...
         4 3 6 4;...
         1 4 6 5];

F_ice = [5 2 8 7;...
         2 3 9 8;...
         3 6 10 9;...
         5 6 10 7];

% Bedrock
patch('parent',H.Ax1,'vertices',V,'faces',F_bed,'facecolor','none','edgecolor',[0.6,0.3,0.0],'linewidth',2);
% Ice
patch('parent',H.Ax1,'vertices',V,'faces',F_ice,'facecolor','none','edgecolor',[0.3,0.5,1.0],'linewidth',2);

set(H.Ax1,'cameraposition',[5.4860   -3.7230    5.5655],'xtick',[],'ytick',[],'ztick',[]);

%% Panel 2: solutions
H.Ax2 = axes('units','pixels','position',[margin_left + wa1 + margin_mid,margin_bot,wa2,ha],'yaxislocation','right',...
  'fontsize',24,'xgrid','on','ygrid','on');

colors = lines(length(foldernames));
% Empty objects for legend
for fi = 1:length(foldernames)
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  line('parent',H.Ax2,'xdata',[],'ydata',[],'color',colors(fi,:),'linewidth',2);
end

% Analytical solution
y = ncread([foldernames{end} '/help_fields_ANT.nc'],'y');
[u_analytical,~] = Schoof2006_SSA_solution( 0.001, 2000, B, L, m, y);

line('parent',H.Ax2,'xdata',y/1e3,'ydata',u_analytical,'linewidth',2,'color','k');
set(H.Ax2,'xlim',[min(y),max(y)]/1e3);
xlabel(H.Ax2,'Distance from centre (km)');
ylabel(H.Ax2,'Ice velocity (m/yr)');

for fi = 1:length(foldernames)
  filename = [foldernames{fi} '/help_fields_ANT.nc'];
  
  time = ncread(filename,'time'); ti = length(time);
  x    = ncread(filename,'x');    xi = round(length(x)/2);
  y    = ncread(filename,'y');
  u    = ncread(filename,'u_base',[xi,1,ti],[1,Inf,1]);
  
  line('parent',H.Ax2,'xdata',y/1e3,'ydata',u,'color',colors(fi,:),'linewidth',2);
end

legend(H.Ax2,'dx = 5 km','dx = 4 km','dx = 3km');
% legend(H.Ax2,'dx = 5 km','dx = 4 km','dx = 3km','dx = 2km','dx = 1km','dx = 500 m','dx = 250 m');