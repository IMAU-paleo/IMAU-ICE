clc
clear all
close all

basefolder = '/Users/berends/Documents/Models/IMAU-ICE/spinup_GRL';

phases = {'phase1_calibration',...
  'phase1a_finecalibration',...
  'phase2_glacialinception',...
  'phase3_termination',...
  'phase4_holocene'};

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
    filename    = [foldername '/scalar_output_global.nc'];
    
    % Safety
    if ~exist( filename,'file')
      % Try the phase1 version
      filename = strrep( filename,'PMIP3ens_','');
      if ~exist( filename,'file')
        continue
      end
    end
    
    results( ri).phase     = phase;
    results( ri).variation = variation;
    results( ri).time      = ncread( filename,'time');
    results( ri).GMSL_GRL  = ncread( filename,'GMSL_GRL');
    
    % Put phase1 in the distant past
    if strcmpi( phase,'phase1_calibration')
      results( ri).time = results( ri).time - results( ri).time( end) - 120000 - 20000;
    elseif strcmpi( phase,'phase1a_finecalibration')
      results( ri).time = results( ri).time - results( ri).time( end) - 120000;
    end
    
  end
end

%% Plot
wa = 900;
ha = 500;

margin_left   = 90;
margin_right  = 25;
margin_bottom = 75;
margin_top    = 50;

wf = margin_left   + wa + margin_right;
hf = margin_bottom + ha + margin_top;

% xlim = [-240,10];
xlim = [-120,0.01];
ylim = [-0.2,1.4];

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margin_left,margin_bottom,wa,ha],'fontsize',24,...
  'xgrid','on','ygrid','on','xlim',xlim,'ylim',ylim);

xlabel(H.Ax,'Time (kyr)')
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
  
  if strcmpi(r.phase,'phase1_calibration') || strcmpi(r.phase,'phase1a_finecalibration')
    linestyle = ':';
  else
    linestyle = '-';
  end
  
  line('parent',H.Ax,'xdata',r.time / 1e3,'ydata',-r.GMSL_GRL,'linewidth',3,'color',colour,'linestyle',linestyle)
  
end

% Mark phases
% line('parent',H.Ax,'xdata',[-120,-120],'ydata',ylim,'linewidth',3,'color','k','linestyle','--')
line('parent',H.Ax,'xdata',[ -21, -21],'ydata',ylim,'linewidth',3,'color','k','linestyle','--')
line('parent',H.Ax,'xdata',[ -10, -10],'ydata',ylim,'linewidth',3,'color','k','linestyle','--')

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