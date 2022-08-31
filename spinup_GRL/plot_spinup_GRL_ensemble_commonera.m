clc
clear all
close all

basefolder = '/Users/berends/Documents/Models/IMAU-ICE/spinup_GRL';

phases = {...
  'phase4_holocene',...
  'phase5_historicalperiod'};

variations = {...
  'hybrid_ZoetIverson_PMIP3ens_40km',...
  'hybrid_ZoetIverson_PMIP3ens_30km',...
  'hybrid_ZoetIverson_PMIP3ens_20km',...
  'hybrid_ZoetIverson_PMIP3ens_16km',...
  'hybrid_ZoetIverson_PMIP3ens_10km',...
  'hybrid_ZoetIverson_HadCM3_20km',...
  'hybrid_ZoetIverson_CCSM_20km',...
  'hybrid_Coulombreg_PMIP3ens_20km'};

%% Read data
results = [];
for phi = 1: length( phases)
  phase = phases{ phi};
  for vi = 1: length( variations)
    
    ri = length( results)+1;
    
    variation   = variations{ vi};
    
    foldername  = [basefolder '/' phase '/' variation];
  
    henk = exist( foldername,'dir');
    if (~(henk== 7)); continue; end
    
    filename    = [foldername '/scalar_output_global.nc'];
    results( ri).phase     = phase;
    results( ri).variation = variation;
    results( ri).time      = ncread( filename,'time');
    results( ri).GMSL_GRL  = ncread( filename,'GMSL_GRL');
    
    % Express time in years CE
%     if strcmpi( phase,'phase4_holocene')
      results( ri).time = results( ri).time + 1950;
%     end
    
  end
end

%% Plot
wa = 900;
ha = 500;

margin_left   = 100;
margin_right  = 25;
margin_bottom = 75;
margin_top    = 50;

wf = margin_left   + wa + margin_right;
hf = margin_bottom + ha + margin_top;

% xlim = [-220,10];
xlim = [0,2050];
ylim = [-0.15,0.15];

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],'fontsize',24,...
  'xgrid','on','ygrid','on','xlim',xlim,'ylim',ylim);

xlabel(H.Ax,'Time (yr CE)')
ylabel(H.Ax,'\Delta VAF (msle)')

variation_colours = linspecer( length( variations),'sequential');

% Empty line objects for legend
legendstr = variations;
for vi = 1: length( variations)
  legendstr{ vi} = parse_underscores( legendstr{ vi});
  colour = variation_colours( vi,:);
  line('parent',H.Ax,'xdata',[],'ydata',[],'linewidth',3,'color',colour)
end

% Plot results
for ri = 1: length( results)
  
  r = results( ri);
  vi = 0;
  for vii = 1: length( variations)
    if strcmpi( variations{ vii},r.variation)
      vi = vii;
      break
    end
  end
  colour = variation_colours( vi,:);
  
  line('parent',H.Ax,'xdata',r.time,'ydata',-r.GMSL_GRL,'linewidth',3,'color',colour)
  
end

% Legend
legend(H.Ax,legendstr,'location','northwest')

function str = parse_underscores( str)
  i = 1;
  while i <= length( str)
    if strcmp( str( i),'_')
      str = [str( 1:i-1) '\_' str( i+1:end)];
      i = i+2;
    else
      i = i+1;
    end
  end
end