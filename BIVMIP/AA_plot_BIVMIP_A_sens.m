clc
clear all
close all

%% Read data
perfect_run = 'BIVMIP_A_perfect_10km';

foldernames = {...
  'BIVMIP_A_inv_10km_perfect',...
  'BIVMIP_A_inv_10km_visc_lo',...
  'BIVMIP_A_inv_10km_visc_hi',...
  'BIVMIP_A_inv_10km_SMB_lo',...
  'BIVMIP_A_inv_10km_SMB_hi',...
  'BIVMIP_A_inv_10km_p_hi',...
  'BIVMIP_A_inv_10km_p_lo',...
  'BIVMIP_A_inv_10km_ut_hi',...
  'BIVMIP_A_inv_10km_ut_lo'};

% Perfect run
filename_restart     = [perfect_run '/restart_ANT.nc'];
filename_help_fields = [perfect_run '/help_fields_ANT.nc'];

x    = ncread(filename_restart,'x');
y    = ncread(filename_restart,'y');
time = ncread(filename_restart,'time');
ti   = length(time);

Hi_perfect  = ncread( filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
phi_perfect = ncread( filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);

% Perturbed runs
results  = [];
dHi_min  =  Inf;
dHi_max  = -Inf;
dphi_min =  Inf;
dphi_max = -Inf;
for fi = 1: length(foldernames)
  
  filename_restart      = [foldernames{fi} '/restart_ANT.nc'];
  filename_help_fields  = [foldernames{fi} '/help_fields_ANT.nc'];
  
  results(fi).time = ncread(filename_restart    ,'time');
  ti = length(results(fi).time);
  
  results(fi).Hi   = ncread(filename_restart    ,'Hi'      ,[1,1,ti],[Inf,Inf,1]);
  results(fi).phi  = ncread(filename_help_fields,'phi_fric',[1,1,ti],[Inf,Inf,1]);
  
  results(fi).dHi  = results(fi).Hi  - Hi_perfect;
  results(fi).dphi = results(fi).phi - phi_perfect;
  
  dHi_min  = min( dHi_min , min(results(fi).dHi( :)));
  dHi_max  = max( dHi_max , max(results(fi).dHi( :)));
  dphi_min = min( dphi_min, min(results(fi).dphi(:)));
  dphi_max = max( dphi_max, max(results(fi).dphi(:)));
end

%% Idea: make histograms of errors in ice thickness / bed roughness

nbins     = 1000;
dHi_bins  = linspace( dHi_min,  dHi_max,  nbins)';
dphi_bins = linspace( dphi_min, dphi_max, nbins)';

for fi = 1: length(foldernames)
  
  % Ice thickness
  results(fi).dHi_hist = zeros( nbins,1);
  biHi = round( (nbins-1) * (results(fi).dHi - dHi_min) / (dHi_max - dHi_min));
  for bi = 1: nbins
    results(fi).dHi_hist( bi) = results(fi).dHi_hist( bi) + sum( biHi(:)==bi);
  end
  
  % Bed roughness
  results(fi).dphi_hist = zeros( nbins,1);
  biphi = round( (nbins-1) * (results(fi).dphi - dphi_min) / (dphi_max - dphi_min));
  for bi = 1: nbins
    results(fi).dphi_hist( bi) = results(fi).dphi_hist( bi) + sum( biphi(:)==bi);
  end
  
  % Apply some smoothing
  results(fi).dHi_hist  = imgaussfilt( results(fi).dHi_hist ,1);
  results(fi).dphi_hist = imgaussfilt( results(fi).dphi_hist,1);
end

%% Plot histograms

wa = 500;
ha = 300;

margins_hor = [75,50,25];
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
    H.Ax( i,j) = axes('parent',H.Fig,'units','pixels','position',[x,y,wa,ha],'fontsize',24);
  end
end

ylabel(H.Ax(1,1),'Relative abundance');
xlabel(H.Ax(1,1),'Error in ice thickness (m)');
xlabel(H.Ax(1,2),['Error in bed roughness (' char(176) ')'])

set(H.Ax(1,1),'ytick',[],'yscale','log');
set(H.Ax(1,2),'ytick',[],'yscale','log');

%% Plot results

% "Perfect" inversion
line('parent',H.Ax(1,1),'xdata',dHi_bins ,'ydata',results(1).dHi_hist ,'color','k','linewidth',2);
line('parent',H.Ax(1,2),'xdata',dphi_bins,'ydata',results(1).dphi_hist,'color','k','linewidth',2);

% "Perturbed" inversions
colors = lines( length(foldernames)-1);
for fi = 2: length(foldernames)
  line('parent',H.Ax(1,1),'xdata',dHi_bins ,'ydata',results(fi).dHi_hist ,'color',colors(fi-1,:),'linewidth',2);
  line('parent',H.Ax(1,2),'xdata',dphi_bins,'ydata',results(fi).dphi_hist,'color',colors(fi-1,:),'linewidth',2);
end

ylim1 = get(H.Ax(1,1),'ylim');
ylim2 = get(H.Ax(1,2),'ylim');
ymax = max(ylim1(2),ylim2(2));
set(H.Ax(1,1),'ylim',[1,ymax]);
set(H.Ax(1,2),'ylim',[1,ymax]);

% Legend
legend(H.Ax(1,1),'Perfect','visc\_lo','visc\_hi','SMB\_lo','SMB\_hi','p\_lo','p\_hi','u_t\_lo','u_t\_hi')