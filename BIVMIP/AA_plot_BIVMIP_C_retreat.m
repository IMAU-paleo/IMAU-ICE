clc
clear all
close all

%% Read data
foldernames = {...
  'BIVMIP_C_retreat_2km_perfect',...
  'BIVMIP_C_retreat_2km_visc_hi',...
  'BIVMIP_C_retreat_2km_visc_lo',...
  'BIVMIP_C_retreat_2km_SMB_hi',...
  'BIVMIP_C_retreat_2km_SMB_lo'};

ice_density      = 910;
seawater_density = 1028;

for fi = 1: length(foldernames)
  
  filename_restart       = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields   = [foldernames{fi} '/help_fields_ANT.nc'];
  
  results(fi).x= ncread(filename_restart    ,'x');
  results(fi).y = ncread(filename_restart    ,'y');
  
  j = find(results(fi).y==0);
  
  results(fi).time = ncread(filename_restart    ,'time');
  results(fi).xGL = zeros( size( results(fi).time));
  
  for ti = 1: length( results(fi).time)
  
    Hi  = ncread( filename_restart,'Hi',[1,j,ti],[Inf,1,1]);
    Hb  = ncread( filename_restart,'Hb',[1,j,ti],[Inf,1,1]);
    TAF = Hi - max(0, (0 - Hb) * (seawater_density / ice_density));
  
    i2 = 2;
    while TAF(i2)>0; i2 = i2+1; end
    i1 = i2-1;
    lambda = TAF(i1) / (TAF(i1) - TAF(i2));
    results(fi).xGL( ti) = lambda * results(fi).x(i2) + (1 - lambda) * results(fi).x(i1);
  end
  
end

%% Plot
wa = 500;
ha = 350;

margins_hor = [100,25];
margins_ver = [80,25];

nax = length(margins_hor)-1;
nay = length(margins_ver)-1;

wf = sum(margins_hor) + nax * wa;
hf = sum(margins_ver) + nay * ha;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = zeros( nay,nax);
H.Axa = zeros( nay,nax);
for i = 1: nay
  for j = 1: nax
    x = sum(margins_hor(1:j )) + (j -1)*wa;
    ip = nay+1-i;
    y = sum(margins_ver(1:ip)) + (ip-1)*ha;
    H.Ax( i,j) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],...
      'fontsize',24,'xgrid','on','ygrid','on');
  end
end

xlabel(H.Ax(1,1),'Time (yr)');
ylabel(H.Ax(1,1),'x_{GL} (km)');

colors = lines(length(foldernames));

for fi = 1: length(foldernames)
  line('parent',H.Ax(1,1),'xdata',results(fi).time,'ydata',(results(fi).xGL/1e3)+400,'linewidth',3,'color',colors(fi,:));
end

legend('Perfect','visc hi','visc lo','SMB hi','SMB lo','location','northeast')