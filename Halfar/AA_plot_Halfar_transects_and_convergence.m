clc
clear all
close all

foldernames = {...
  '/Users/berends/Documents/Models/IMAU-ICE/Halfar/Halfar_50km',...
  '/Users/berends/Documents/Models/IMAU-ICE/Halfar/Halfar_25km',...
  '/Users/berends/Documents/Models/IMAU-ICE/Halfar/Halfar_16km',...
  '/Users/berends/Documents/Models/IMAU-ICE/Halfar/Halfar_8km'};

linestyles = {'-','--',':','-.'};

for fi = 1: length(foldernames)
  filename = [foldernames{fi} '/restart_ANT.nc'];
  results(fi).x  = ncread(filename,'x');
  y  = ncread(filename,'y');
  yi = find(y==0);
  results(fi).Hi   = ncread(filename,'Hi',[1,yi,1],[Inf,1,Inf]);
  results(fi).Time = ncread(filename,'time');
end

times_to_plot = [0,10,100]*1e3;

%% Ice thickness

wa = 400;
ha = 300;

margin_left   = 100;
margin_middle = 25;
margin_right  = 25;
margin_bottom = 80;
margin_top    = 50;

wf = margin_left + wa + margin_middle + wa + margin_middle + wa + margin_right;
hf = margin_bottom + ha + margin_top;

H1.Fig = figure('color','w','position',[400,300,wf,hf]);
H1.Ax(1) = axes('parent',H1.Fig,'units','pixels','position',[margin_left+(wa+margin_middle)*0,margin_bottom,wa,ha]);
H1.Ax(2) = axes('parent',H1.Fig,'units','pixels','position',[margin_left+(wa+margin_middle)*1,margin_bottom,wa,ha]);
H1.Ax(3) = axes('parent',H1.Fig,'units','pixels','position',[margin_left+(wa+margin_middle)*2,margin_bottom,wa,ha]);

for a = 1:3
  set(H1.Ax(a),'fontsize',24,'xgrid','on','ygrid','on','xlim',[0,750],'ylim',[0,5000],'xtick',0:100:700,'ytick',0:500:5000)
end
set(H1.Ax(1),'yticklabels',{'0','','1000','','2000','','3000','','4000','','5000'});
set(H1.Ax(2),'yticklabels',[]);
set(H1.Ax(3),'yticklabels',[]);

xlabel(H1.Ax(2),'Distance from divide (km)')
ylabel(H1.Ax(1),'Ice thickness (m)')
title(H1.Ax(1),'a')
title(H1.Ax(2),'b')
title(H1.Ax(3),'c')

for a=1:3
  text(H1.Ax(a),530,4750,['t = ' num2str(times_to_plot(a)/1e3) ' kyr'],'fontsize',24,'fontweight','bold');
end

for a = 1:3

  Halfar.x = linspace(0,max(results(1).x),10000);
  Halfar.Hi = Halfar_solution(5000,300000,Halfar.x,0,times_to_plot(a));
  line('parent',H1.Ax(a),'xdata',Halfar.x/1e3,'ydata',Halfar.Hi,'color','k','linewidth',1);

  for fi = 1: length(foldernames)
    ti = find(results(fi).Time == times_to_plot(a));
    Hi = results(fi).Hi(:,1,ti);
    line('parent',H1.Ax(a),'xdata',results(fi).x/1e3,'ydata',Hi,'color','b','linewidth',2,'linestyle',linestyles{fi});
  end
end

%% Convergence in margin position

wa = 500;
ha = 500;
margin_left = 100;
margin_bottom = 100;
margin_right = 25;
margin_top = 25;
wf = margin_left + wa + margin_right;
hf = margin_bottom + ha + margin_top;

H2.Fig = figure('color','w','position',[300,300,wf,hf]);
H2.Ax  = axes('parent',H2.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],...
  'xscale','log','yscale','log','fontsize',24,'xgrid','on','ygrid','on');
xlabel(H2.Ax,'Resolution (km)')
ylabel(H2.Ax,'Margin error (km)')

% Empty line object for legend
line('parent',H2.Ax,'xdata',[],'ydata',[],'linestyle','--','color','b','linewidth',2);

resolutions = [50,25,16,8,4];
margin_error = zeros(length(foldernames),1);

x_exact = linspace(0,max(results(1).x),100000);
Hi_exact = Halfar_solution(5000,300000,x_exact,0,times_to_plot(3));
xi = length(x_exact);
while (Hi_exact(xi)==0); xi=xi-1; end
margin_exact = x_exact(xi);

for fi = 1: length(foldernames)
  
  ti = find(results(fi).Time == times_to_plot(3));
  Hi = results(fi).Hi(:,1,ti);
  xi = length(results(fi).x);
  while (Hi(xi)==0); xi=xi-1; end
  margin_mod = results(fi).x(xi);
  
  margin_error(fi) = (margin_mod - margin_exact)/1e3;
  
end

line('parent',H2.Ax,'xdata',resolutions(1:length(foldernames)),'ydata',margin_error,'linestyle','none',...
  'marker','o','markerfacecolor','b','markeredgecolor','b','markersize',10);

% Loglinear fit
p = polyfit( log(resolutions(1:length(foldernames))), log(margin_error(1:length(foldernames))), 1);
RMSE_fit = exp(polyval(p,log(resolutions(1:length(foldernames)))));
line('parent',H2.Ax,'xdata',resolutions(1:length(foldernames)),'ydata',RMSE_fit,'color','b','linestyle','--','linewidth',2)

legend(H2.Ax,['O(R^{' num2str(round(p(1)*100)/100) '})'],'location','northwest');