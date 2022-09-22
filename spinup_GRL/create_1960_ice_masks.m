clc
clear all
close all

icemasks_foldername = 'IMAUICE_icemasks_1960';

foldernames = {...
  'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_40km',...
  'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_30km',...
  'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_20km',...
  'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_16km',...
  'phase4_holocene/hybrid_ZoetIverson_PMIP3ens_10km',...
  'phase4_holocene/hybrid_ZoetIverson_HadCM3_20km',...
  'phase4_holocene/hybrid_ZoetIverson_CCSM_20km',...
  'phase4_holocene/hybrid_Coulombreg_PMIP3ens_20km'};

if exist(icemasks_foldername,'dir')
  henk = dir(icemasks_foldername);
  for i = 3:length(henk)
    delete([icemasks_foldername '/' henk(i).name])
  end
  rmdir(icemasks_foldername)
end
mkdir(icemasks_foldername)

fid = fopen([icemasks_foldername '/readme.txt'],'w');

for fi = 1: length( foldernames)
  
  foldername = foldernames{ fi};

  filename_src = [foldername '/restart_GRL.nc'];
  filename_dst = [icemasks_foldername '/IMAUICE' num2str(fi) '_sftgif.nc'];

  %% Read data

  grid.x = ncread( filename_src,'x'); grid.nx = length( grid.x);
  grid.y = ncread( filename_src,'y'); grid.ny = length( grid.y);

  time   = ncread( filename_src,'time'); ti = find( time==10);
  ice.Hi = ncread( filename_src,'Hi',[1,1,ti],[Inf,Inf,1]);
  sftgif = zeros( size( ice.Hi));
  sftgif( ice.Hi > 0.1) = 1;

  %% Set up NetCDF template

  % General
  f.Filename         = filename_dst;
  f.Name             = '/';
  f.Groups           = [];
  f.Format           = 'classic';
  f.Attributes       = [];

  % Dimensions
  dim_x.Name         = 'x';
  dim_x.Length       = grid.nx;
  dim_x.Unlimited    = false;

  dim_y.Name         = 'y';
  dim_y.Length       = grid.ny;
  dim_y.Unlimited    = false;

  f.Dimensions(1)    = dim_x;
  f.Dimensions(2)    = dim_y;

  % x,y variables
  var_x.Name         = 'x';
  var_x.Dimensions   = dim_x;
  var_x.Size         = dim_x.Length;
  var_x.Datatype     = 'double';
  var_x.Attributes(1).Name  = 'long_name';
  var_x.Attributes(1).Value = 'x-coordinate';
  var_x.Attributes(2).Name  = 'units';
  var_x.Attributes(2).Value = 'm';
  var_x.ChunkSize    = [];
  var_x.DeflateLevel = [];
  var_x.Shuffle      = false;

  var_y = var_x;
  var_y.Name         = 'y';
  var_y.Dimensions   = dim_y;
  var_y.Size         = dim_y.Length;
  var_y.Attributes(1).Value = 'y-coordinate';

  f.Variables(1) = var_x;
  f.Variables(2) = var_y;

  % field variable
  var.Name          = 'sftgif';
  var.Dimensions(1) = dim_x;
  var.Dimensions(2) = dim_y;
  var.Size          = [dim_x.Length, dim_y.Length];
  var.Datatype      = 'double';
  var.Attributes(1).Name  = 'long_name';
  var.Attributes(1).Value = 'land_ice_area_fraction';
  var.Attributes(2).Name  = 'units';
  var.Attributes(2).Value = '1';
  var.ChunkSize     = [];
  var.DeflateLevel  = [];
  var.Shuffle       = false;

  f.Variables(3) = var;

  % Create empty NetCDF file
  if exist( filename_dst,'file')
    delete( filename_dst)
  end

  ncwriteschema( filename_dst, f);

  % Write data
  ncwrite( filename_dst,'x'     ,grid.x);
  ncwrite( filename_dst,'y'     ,grid.y);
  ncwrite( filename_dst,'sftgif',sftgif);
  
  % Save name of run in readme
  spinupname = strrep(foldername,'phase4_holocene/','');
  fprintf(fid,['file "' filename_dst '" corresponds to spin-up variation "' spinupname '"\n']);

end

fclose(fid);