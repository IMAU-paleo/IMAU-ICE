clc
clear all
close all

foldernames = {...
  'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_40km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_30km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_20km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_16km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_10km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_HadCM3_20km',...
  'phase5_historicalperiod/hybrid_ZoetIverson_CCSM_20km',...
  'phase5_historicalperiod/hybrid_Coulombreg_PMIP3ens_20km'};

tstart = 10;
tstop  = 39;

for fi = 1: length( foldernames)
  foldername = foldernames{ fi};

  %% Calculate average temperature and SMB over this period
  filename = [foldername '/help_fields_GRL.nc'];
  time = ncread( filename, 'time');
  ti_start = 1;
  while time( ti_start) < tstart
    ti_start = ti_start + 1;
  end
  ti_stop = length( time);
  while time( ti_stop) > tstop
    ti_stop = ti_stop - 1;
  end

  x   = ncread( filename,'x'); nx = length( x);
  y   = ncread( filename,'y'); ny = length( y);
  T2m = ncread( filename,'T2m_year',[1,1,ti_start],[Inf,Inf,ti_stop+1-ti_start]);
  SMB = ncread( filename,'SMB_year',[1,1,ti_start],[Inf,Inf,ti_stop+1-ti_start]);
  Hs  = ncread( filename,'Hs'      ,[1,1,ti_start],[Inf,Inf,ti_stop+1-ti_start]);

  % Average SMB and Hs
  SMB_av = mean( SMB,3);
  Hs_av  = mean( Hs ,3);

  % Apply a simple lapse-rate correction to account for changes in ice geometry
  lambda = 0.008;
  Hs  = ncread( filename,'Hs',[1,1,ti_start],[Inf,Inf,ti_stop+1-ti_start]);
  T2m_SL = T2m + Hs * lambda;
  T2m_SL_av = mean( T2m_SL,3);
  T2m_av = T2m_SL_av - Hs(:,:,end) * lambda;

  %% Set up a NetCDF template
  
  filename_dst = [foldername '/baseline_climate_1960_1989.nc'];

  % General
  f.Filename             = filename_dst;
  f.Name                 = '/';
  f.Attributes           = [];
  f.Groups               = [];
  f.Format               = 'classic';

  % Dimensions
  dim_x.Name             = 'x';
  dim_x.Length           = nx;
  dim_x.Unlimited        = false;

  dim_y.Name             = 'y';
  dim_y.Length           = ny;
  dim_y.Unlimited        = false;

  f.Dimensions( 1) = dim_x;
  f.Dimensions( 2) = dim_y;

  % Variables
  var_x.Name             = 'x';
  var_x.Dimensions       = dim_x;
  var_x.Size             = nx;
  var_x.Datatype         = 'double';
  var_x.ChunkSize        = [];
  var_x.FillValue        = [];
  var_x.DeflateLevel     = [];
  var_x.Shuffle          = false;
  var_x.Attributes( 1).Name  = 'long_name';
  var_x.Attributes( 1).Value = 'x-coordinate';
  var_x.Attributes( 2).Name  = 'units';
  var_x.Attributes( 2).Value = 'm';

  var_y.Name             = 'y';
  var_y.Dimensions       = dim_y;
  var_y.Size             = ny;
  var_y.Datatype         = 'double';
  var_y.ChunkSize        = [];
  var_y.FillValue        = [];
  var_y.DeflateLevel     = [];
  var_y.Shuffle          = false;
  var_y.Attributes( 1).Name  = 'long_name';
  var_y.Attributes( 1).Value = 'y-coordinate';
  var_y.Attributes( 2).Name  = 'units';
  var_y.Attributes( 2).Value = 'm';

  var_SMB.Name             = 'SMB';
  var_SMB.Dimensions( 1)   = dim_x;
  var_SMB.Dimensions( 2)   = dim_y;
  var_SMB.Size             = [nx,ny];
  var_SMB.Datatype         = 'double';
  var_SMB.ChunkSize        = [];
  var_SMB.FillValue        = [];
  var_SMB.DeflateLevel     = [];
  var_SMB.Shuffle          = false;
  var_SMB.Attributes( 1).Name  = 'long_name';
  var_SMB.Attributes( 1).Value = 'Surface mass balance';
  var_SMB.Attributes( 2).Name  = 'units';
  var_SMB.Attributes( 2).Value = 'm.i.e./yr';

  var_ST.Name             = 'ST';
  var_ST.Dimensions( 1)   = dim_x;
  var_ST.Dimensions( 2)   = dim_y;
  var_ST.Size             = [nx,ny];
  var_ST.Datatype         = 'double';
  var_ST.ChunkSize        = [];
  var_ST.FillValue        = [];
  var_ST.DeflateLevel     = [];
  var_ST.Shuffle          = false;
  var_ST.Attributes( 1).Name  = 'long_name';
  var_ST.Attributes( 1).Value = 'Annual mean 2-m air temperature';
  var_ST.Attributes( 2).Name  = 'units';
  var_ST.Attributes( 2).Value = 'K';

  var_Hs.Name             = 'Hs';
  var_Hs.Dimensions( 1)   = dim_x;
  var_Hs.Dimensions( 2)   = dim_y;
  var_Hs.Size             = [nx,ny];
  var_Hs.Datatype         = 'double';
  var_Hs.ChunkSize        = [];
  var_Hs.FillValue        = [];
  var_Hs.DeflateLevel     = [];
  var_Hs.Shuffle          = false;
  var_Hs.Attributes( 1).Name  = 'long_name';
  var_Hs.Attributes( 1).Value = 'Surface elevation';
  var_Hs.Attributes( 2).Name  = 'units';
  var_Hs.Attributes( 2).Value = 'm';

  f.Variables( 1) = var_x;
  f.Variables( 2) = var_y;
  f.Variables( 3) = var_SMB;
  f.Variables( 4) = var_ST;
  f.Variables( 5) = var_Hs;

  %% Create NetCDF file

  if exist( filename_dst,'file')
    error(['File "' filename_dst '" already exists!'])
  end

  ncwriteschema( filename_dst, f);

  ncwrite( filename_dst, 'x'  , x     );
  ncwrite( filename_dst, 'y'  , y     );
  ncwrite( filename_dst, 'ST' , T2m_av);
  ncwrite( filename_dst, 'SMB', SMB_av);
  ncwrite( filename_dst, 'Hs' , Hs_av );

end
