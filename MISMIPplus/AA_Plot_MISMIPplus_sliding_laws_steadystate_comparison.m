clc
clear all
close all

foldernames = {...
  'MISMIPplus_init_5km_slid1_DIVA',...
  'MISMIPplus_init_5km_slid2_DIVA',...
  'MISMIPplus_init_5km_slid3_DIVA',...
  'MISMIPplus_init_5km_slid1_SIASSA',...
  'MISMIPplus_init_5km_slid2_SIASSA',...
  'MISMIPplus_init_5km_slid3_SIASSA',...
  'MISMIPplus_init_2km_slid1_DIVA',...
  'MISMIPplus_init_2km_slid2_DIVA',...
  'MISMIPplus_init_2km_slid3_DIVA'};

%% Read model output
for fi = 1:length(foldernames)
  
  filename = [foldernames{fi} '/restart_ANT.nc'];
  time = ncread(filename,'time');
  ti = length(time);
  
  x = ncread(filename,'x');
  y = ncread(filename,'y');
  j = ceil(length(y)/2);
  
  results(fi).Hi  = ncread(filename,'Hi',[1,j,ti],[Inf,1,1]);
  results(fi).Hb  = ncread(filename,'Hb',[1,j,ti],[Inf,1,1]);
  results(fi).Hs  = results(fi).Hi + max( -910 / 1028 * results(fi).Hi, results(fi).Hb);
  
  results(fi).Hs( results(fi).Hi==0) = NaN;
  results(fi).Hi( results(fi).Hi==0) = NaN;
  
  i0 = length(x); while isnan(results(fi).Hi(i0,1)); i0 = i0-1; end
  results(fi).Hi(i0+1,:) = 0;
  results(fi).Hs(i0+1,:) = 0;
  
  results(fi).x   = (x - min(x))/range(x) * 800;
  results(fi).Hib = results(fi).Hs - results(fi).Hi;
  
end

%% Plot

wa = 600;
ha = 400;

margin_left   = 110;
margin_right  = 25;
margin_bottom = 80;
margin_top    = 25;

wf = margin_left   + wa + margin_right;
hf = margin_bottom + ha + margin_top;

H.Fig = figure('position',[300,300,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],'fontsize',24,...
  'xgrid','on','ygrid','on','xlim',[0,800],'ylim',[-800,2000]);
xlabel(H.Ax,'x (km)')
ylabel(H.Ax,'z (m)')

% empty line objects for legend
colors = lines(length(foldernames));

for fi = 1:length(foldernames)
  line('xdata',results(fi).x,'ydata',results(fi).Hs, 'linewidth',2,'color',colors(fi,:));
  line('xdata',results(fi).x,'ydata',results(fi).Hib,'linewidth',2,'color',colors(fi,:));
end
line('xdata',results(end).x,'ydata',results(end).Hb, 'linestyle','-','linewidth',2,'color','k');

% legend(H.Ax,'Sliding law #1','Sliding law #2','Sliding law #3','location','northeast');