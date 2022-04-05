clc
clear all
close all

% Create NetCDF files containing the bed roughness, bed topography, and SMB
% (including perturbed versions) for experiment I (based on the EISMINT-I
% "moving margin" experiment, but slightly smaller and with sliding).

%% Define parameters

% The resolutions at which we want to have the data
resolutions = [40,20,10] * 1e3;

% Domain size
xmin        = -700*1e3;     % x-coordinate of western  domain border                    [m]
xmax        =  700*1e3;     % x-coordinate of eastern  domain border                    [m]
ymin        = -700*1e3;     % y-coordinate of southern domain border                    [m]
ymax        =  700*1e3;     % y-coordinate of northern domain border                    [m]

% Bed roughness parameters
phi_min     = 0.1;          % Till friction angle in the centre of the ice stream       [degrees]
phi_max     = 5.0;          % Till friction angle outside of the ice stream             [degrees]
x_c         = 0;            % x-coordinate of ice-stream centre                         [m]
y_c         = -400*1e3;     % y-coordinate of ice-stream centre                         [m]
sigma_x     =   50*1e3;     % x-direction ice-stream half-width                         [m]
sigma_y     =  300*1e3;     % y-direction ice-stream half-width                         [m]

% Topography parameters
A_topo_hi   =  250;         % High-perturbed bed elevation error at the ice divide      [m]
A_topo_lo   = -250;         % Low -perturbed bed elevation error at the ice divide      [m]
sigma_topo  = 450*1e3;      % Radius of bed elevation error                             [m]

% SMB parameters
M_max       = 0.5;          % Maximum accumulation rate                                 [m/yr]
E           = 400*1e3;      % Radius of accumulation zone                               [m]
S           = 1e5;          % Melt rate increase over radial distance from grid centre  [yr]
f_SMB_hi    = 1.1;          % Multiplication factor for high perturbed SMB              [-]
f_SMB_lo    = 0.9;          % Multiplication factor for low  perturbed SMB              [-]

% File names
filename_bed_roughness  = 'exp_I_bed_roughness';  % Will be extended with the resolution and file type, e.g. "exp_I_bed_roughness_40km.nc"
filename_bed_topography = 'exp_I_bed_topography'; % In the case of topography and SMB, an additional extension for the two perturbed versions
filename_SMB            = 'exp_I_SMB';            % will be added, e.g. "exp_I_SMB_40km_hi.nc".

for ri = 1:length( resolutions)
  
  %% Define the grid
  
  % The grid resolution
  resolution = resolutions( ri);
  
  % The resolution filename extension
  str_res = ['_' num2str( resolution / 1e3) 'km'];
  
  % Create the grid
  grid = create_grid( xmin, xmax, ymin, ymax, resolution);
  
  %% Calculate data
  
  % Bed roughness
  phi = zeros( grid.nx, grid.ny);
  
  for i = 1: grid.nx
    for j = 1: grid.ny
      phi( i,j) = calc_bed_roughness( phi_min, phi_max, x_c, y_c, sigma_x, sigma_y, grid.x( i), grid.y( j));
    end
  end
  
  % Bed topography
  b    = zeros( grid.nx, grid.ny);
  b_hi = zeros( grid.nx, grid.ny);
  b_lo = zeros( grid.nx, grid.ny);
  
  for i = 1: grid.nx
    for j = 1: grid.ny
      
      % Regular version: just a flat earth
      b( i,j) = 0;
      
      % Perturbed versions
      db_hi = calc_topography_error_exp_I( A_topo_hi, sigma_topo, grid.x( i), grid.y( j));
      db_lo = calc_topography_error_exp_I( A_topo_lo, sigma_topo, grid.x( i), grid.y( j));
      
      b_hi( i,j) = b( i,j) + db_hi;
      b_lo( i,j) = b( i,j) + db_lo;
      
    end
  end
  
  % Surface mass balance
  SMB    = zeros( grid.nx, grid.ny);
  SMB_hi = zeros( grid.nx, grid.ny);
  SMB_lo = zeros( grid.nx, grid.ny);
  
  for i = 1: grid.nx
    for j = 1: grid.ny
      
      % Regular version
      SMB( i,j) = calc_SMB_exp_I( M_max, E, S, grid.x( i), grid.y( j));
      
      % Perturbed versions
      SMB_hi( i,j) = SMB( i,j) * f_SMB_hi;
      SMB_lo( i,j) = SMB( i,j) * f_SMB_lo;
      
    end
  end
  
  %% Create and write to NetCDF files
  
  % =========================
  % ===== Bed roughness =====
  % =========================
  
  filename = [filename_bed_roughness str_res '.nc'];
  
  % Delete existing file
  if exist( filename,'file')
    delete( filename)
  end
  
  % NetCDF template
  f = create_NetCDF_template_bed_roughness( grid, filename);
  
  % Create file
  ncwriteschema( filename, f);
  
  % Write data
  ncwrite( filename,'x'       ,grid.x);
  ncwrite( filename,'y'       ,grid.y);
  ncwrite( filename,'phi_fric',phi   );
  
  % ==========================
  % ===== Bed topography =====
  % ==========================
  
  filename    = [filename_bed_topography str_res '.nc'];
  filename_hi = [filename_bed_topography str_res '_hi.nc'];
  filename_lo = [filename_bed_topography str_res '_lo.nc'];
  
  % Delete existing files
  if exist( filename,'file')
    delete( filename)
  end
  if exist( filename_hi,'file')
    delete( filename_hi)
  end
  if exist( filename_lo,'file')
    delete( filename_lo)
  end
  
  % NetCDF templates
  f = create_NetCDF_template_bed_topography( grid, filename);
  
  f_hi = f;
  f_lo = f;
  f_hi.Filename = filename_hi;
  f_lo.Filename = filename_lo;
  
  % Create files
  ncwriteschema( filename   , f   );
  ncwriteschema( filename_hi, f_hi);
  ncwriteschema( filename_lo, f_lo);
  
  % Write data
  ncwrite( filename   ,'x'       ,grid.x);
  ncwrite( filename   ,'y'       ,grid.y);
  ncwrite( filename   ,'Hi'      ,zeros( grid.nx, grid.ny));
  ncwrite( filename   ,'Hb'      ,b     );
  ncwrite( filename   ,'Hs'      ,zeros( grid.nx, grid.ny));
  
  ncwrite( filename_hi,'x'       ,grid.x);
  ncwrite( filename_hi,'y'       ,grid.y);
  ncwrite( filename_hi,'Hi'      ,zeros( grid.nx, grid.ny));
  ncwrite( filename_hi,'Hb'      ,b     );
  ncwrite( filename_hi,'Hs'      ,zeros( grid.nx, grid.ny));
  
  ncwrite( filename_lo,'x'       ,grid.x);
  ncwrite( filename_lo,'y'       ,grid.y);
  ncwrite( filename_lo,'Hi'      ,zeros( grid.nx, grid.ny));
  ncwrite( filename_lo,'Hb'      ,b     );
  ncwrite( filename_lo,'Hs'      ,zeros( grid.nx, grid.ny));
  
  % ================================
  % ===== Surface mass balance =====
  % ================================
  
  filename    = [filename_SMB str_res '.nc'];
  filename_hi = [filename_SMB str_res '_hi.nc'];
  filename_lo = [filename_SMB str_res '_lo.nc'];
  
  % Delete existing files
  if exist( filename,'file')
    delete( filename)
  end
  if exist( filename_hi,'file')
    delete( filename_hi)
  end
  if exist( filename_lo,'file')
    delete( filename_lo)
  end
  
  % NetCDF templates
  f = create_NetCDF_template_SMB( grid, filename);
  
  f_hi = f;
  f_lo = f;
  f_hi.Filename = filename_hi;
  f_lo.Filename = filename_lo;
  
  % Create files
  ncwriteschema( filename   , f   );
  ncwriteschema( filename_hi, f_hi);
  ncwriteschema( filename_lo, f_lo);
  
  % For IMAU-ICE, an SMB forcing file needs a time dimension
  time = [-1e6,1e6];
  SMB_time    = zeros( grid.nx, grid.ny, 2);
  SMB_hi_time = zeros( grid.nx, grid.ny, 2);
  SMB_lo_time = zeros( grid.nx, grid.ny, 2);
  for ti = 1:2
    SMB_time(    :,:,ti) = SMB;
    SMB_hi_time( :,:,ti) = SMB_hi;
    SMB_lo_time( :,:,ti) = SMB_lo;
  end
  
  % It also needs temperature, even though in these experiments it's not
  % actually used...
  
  T2m = zeros( grid.nx, grid.ny, 2);
  
  % Write data
  ncwrite( filename   ,'x'       ,grid.x     );
  ncwrite( filename   ,'y'       ,grid.y     );
  ncwrite( filename   ,'time'    ,time       );
  ncwrite( filename   ,'SMB'     ,SMB_time   );
  ncwrite( filename   ,'T2m'     ,T2m        );
  
  ncwrite( filename_hi,'x'       ,grid.x     );
  ncwrite( filename_hi,'y'       ,grid.y     );
  ncwrite( filename_hi,'time'    ,time       );
  ncwrite( filename_hi,'SMB'     ,SMB_hi_time);
  ncwrite( filename_hi,'T2m'     ,T2m        );
  
  ncwrite( filename_lo,'x'       ,grid.x     );
  ncwrite( filename_lo,'y'       ,grid.y     );
  ncwrite( filename_lo,'time'    ,time       );
  ncwrite( filename_lo,'SMB'     ,SMB_lo_time);
  ncwrite( filename_lo,'T2m'     ,T2m        );
  
end

function grid = create_grid( xmin, xmax, ymin, ymax, dx)
  % Create a square grid
  %
  % Code copied from IMAU-ICE
  
  grid.dx = dx;
      
  % Determine the number of grid cells that can fit in this domain
  xmid = (xmax + xmin) / 2;
  ymid = (ymax + ymin) / 2;
  nsx  = floor( (xmax - xmid) / grid.dx);
  nsy  = floor( (ymax - ymid) / grid.dx);

  grid.nx = 1 + 2*nsx;
  grid.ny = 1 + 2*nsy;
  
  % Fill in x and y
  grid.x = zeros( grid.nx,1);
  grid.y = zeros( grid.ny,1);
  
  for i = 1: grid.nx
    grid.x( i) = -nsx*grid.dx + (i-1)*grid.dx;
  end
  for j = 1: grid.ny
    grid.y( j) = -nsy*grid.dx + (j-1)*grid.dx;
  end
end
function phi  = calc_bed_roughness( phi_min, phi_max, x_c, y_c, sigma_x, sigma_y, x, y)
  % Calculate the spatially variable bed roughness field
  %
  % Basically a uniform value of phi_max, with a single ice stream at [x_c,y_c]
  % where it decreases to phi_min at the centre, with a half-width of
  % sigma_x, sigma_y in the respective x- and y-directions.
  
  phi = phi_max - (phi_max - phi_min) * exp( -0.5 * (((x - x_c) / sigma_x)^2 + ((y - y_c) / sigma_y)^2));
  
end
function db   = calc_topography_error_exp_I( A, sigma, x, y)
  % Calculate the introduced bedrock elevation error
  % 
  % A single sinus-shaped bump (so that the error is zero outside the ice margin)

  r = sqrt( x^2 + y^2);

  if r > sigma
    db = 0;
  else
    db = A * 0.5 * (1 + cos( r * pi / sigma));
  end

end
function SMB  = calc_SMB_exp_I( M_max, E, S, x, y)
  % Calculate the surface mass balance
  %
  % Based on the EISMINT-I "moving margin" experiment, but with a slightly
  % smaller radius
  
  r = sqrt( x^2 + y^2);
  
  SMB = min( M_max, (E - r) / S);
  
end
function f    = create_NetCDF_template_bed_roughness(  grid, filename)
  % Create a template for the bed roughness NetCDF file
  
  % Metadata
  f.Filename   = filename;
  f.Name       = '/';
  
  % Dimensions
  f.Dimensions(1).Name      = 'x';
  f.Dimensions(1).Length    = grid.nx;
  f.Dimensions(1).Unlimited = false;
  
  f.Dimensions(2).Name      = 'y';
  f.Dimensions(2).Length    = grid.ny;
  f.Dimensions(2).Unlimited = false;
  
  % Dimension variables
  
  % x
  f.Variables(1).Name         = 'x';
  f.Variables(1).Dimensions   = f.Dimensions(1);
  f.Variables(1).Size         = grid.nx;
  f.Variables(1).Datatype     = 'double';
  f.Variables(1).Attributes(1).Name  = 'long_name';
  f.Variables(1).Attributes(1).Value = 'X-coordinate';
  f.Variables(1).Attributes(2).Name  = 'units';
  f.Variables(1).Attributes(2).Value = 'm';
  f.Variables(1).ChunkSize    = [];
  f.Variables(1).FillValue    = [];
  f.Variables(1).DeflateLevel = [];
  f.Variables(1).Shuffle      = false;
  
  % y
  f.Variables(2).Name         = 'y';
  f.Variables(2).Dimensions   = f.Dimensions(2);
  f.Variables(2).Size         = grid.ny;
  f.Variables(2).Datatype     = 'double';
  f.Variables(2).Attributes(1).Name  = 'long_name';
  f.Variables(2).Attributes(1).Value = 'Y-coordinate';
  f.Variables(2).Attributes(2).Name  = 'units';
  f.Variables(2).Attributes(2).Value = 'm';
  f.Variables(2).ChunkSize    = [];
  f.Variables(2).FillValue    = [];
  f.Variables(2).DeflateLevel = [];
  f.Variables(2).Shuffle      = false;
  
  % Bed roughness variable
  
  f.Variables(3).Name         = 'phi_fric';
  f.Variables(3).Dimensions   = f.Dimensions;
  f.Variables(3).Size         = [grid.nx, grid.ny];
  f.Variables(3).Datatype     = 'double';
  f.Variables(3).Attributes(1).Name  = 'long_name';
  f.Variables(3).Attributes(1).Value = 'Till friction angle';
  f.Variables(3).Attributes(2).Name  = 'units';
  f.Variables(3).Attributes(2).Value = 'degrees';
  f.Variables(3).ChunkSize    = [];
  f.Variables(3).FillValue    = [];
  f.Variables(3).DeflateLevel = [];
  f.Variables(3).Shuffle      = false;
  
  % Final metadata
  f.Attributes = [];
  f.Groups     = [];
  f.Format     = 'classic';
  
end
function f    = create_NetCDF_template_bed_topography( grid, filename)
  % Create a template for the bed topography NetCDF file
  
  % Metadata
  f.Filename   = filename;
  f.Name       = '/';
  
  % Dimensions
  f.Dimensions(1).Name      = 'x';
  f.Dimensions(1).Length    = grid.nx;
  f.Dimensions(1).Unlimited = false;
  
  f.Dimensions(2).Name      = 'y';
  f.Dimensions(2).Length    = grid.ny;
  f.Dimensions(2).Unlimited = false;
  
  % Dimension variables
  
  % x
  f.Variables(1).Name         = 'x';
  f.Variables(1).Dimensions   = f.Dimensions(1);
  f.Variables(1).Size         = grid.nx;
  f.Variables(1).Datatype     = 'double';
  f.Variables(1).Attributes(1).Name  = 'long_name';
  f.Variables(1).Attributes(1).Value = 'X-coordinate';
  f.Variables(1).Attributes(2).Name  = 'units';
  f.Variables(1).Attributes(2).Value = 'm';
  f.Variables(1).ChunkSize    = [];
  f.Variables(1).FillValue    = [];
  f.Variables(1).DeflateLevel = [];
  f.Variables(1).Shuffle      = false;
  
  % y
  f.Variables(2).Name         = 'y';
  f.Variables(2).Dimensions   = f.Dimensions(2);
  f.Variables(2).Size         = grid.ny;
  f.Variables(2).Datatype     = 'double';
  f.Variables(2).Attributes(1).Name  = 'long_name';
  f.Variables(2).Attributes(1).Value = 'Y-coordinate';
  f.Variables(2).Attributes(2).Name  = 'units';
  f.Variables(2).Attributes(2).Value = 'm';
  f.Variables(2).ChunkSize    = [];
  f.Variables(2).FillValue    = [];
  f.Variables(2).DeflateLevel = [];
  f.Variables(2).Shuffle      = false;
  
  % Ice thickness variable
  
  f.Variables(3).Name         = 'Hi';
  f.Variables(3).Dimensions   = f.Dimensions;
  f.Variables(3).Size         = [grid.nx, grid.ny];
  f.Variables(3).Datatype     = 'double';
  f.Variables(3).Attributes(1).Name  = 'long_name';
  f.Variables(3).Attributes(1).Value = 'Ice thickness';
  f.Variables(3).Attributes(2).Name  = 'units';
  f.Variables(3).Attributes(2).Value = 'm';
  f.Variables(3).ChunkSize    = [];
  f.Variables(3).FillValue    = [];
  f.Variables(3).DeflateLevel = [];
  f.Variables(3).Shuffle      = false;
  
  % Bed topography variable
  
  f.Variables(4).Name         = 'Hb';
  f.Variables(4).Dimensions   = f.Dimensions;
  f.Variables(4).Size         = [grid.nx, grid.ny];
  f.Variables(4).Datatype     = 'double';
  f.Variables(4).Attributes(1).Name  = 'long_name';
  f.Variables(4).Attributes(1).Value = 'Bedrock elevation';
  f.Variables(4).Attributes(2).Name  = 'units';
  f.Variables(4).Attributes(2).Value = 'm';
  f.Variables(4).ChunkSize    = [];
  f.Variables(4).FillValue    = [];
  f.Variables(4).DeflateLevel = [];
  f.Variables(4).Shuffle      = false;
  
  % Surface elevation variable
  
  f.Variables(5).Name         = 'Hs';
  f.Variables(5).Dimensions   = f.Dimensions;
  f.Variables(5).Size         = [grid.nx, grid.ny];
  f.Variables(5).Datatype     = 'double';
  f.Variables(5).Attributes(1).Name  = 'long_name';
  f.Variables(5).Attributes(1).Value = 'Surface elevation';
  f.Variables(5).Attributes(2).Name  = 'units';
  f.Variables(5).Attributes(2).Value = 'm';
  f.Variables(5).ChunkSize    = [];
  f.Variables(5).FillValue    = [];
  f.Variables(5).DeflateLevel = [];
  f.Variables(5).Shuffle      = false;
  
  % Final metadata
  f.Attributes = [];
  f.Groups     = [];
  f.Format     = 'classic';
  
end
function f    = create_NetCDF_template_SMB(            grid, filename)
  % Create a template for the SMB NetCDF file
  
  % Metadata
  f.Filename   = filename;
  f.Name       = '/';
  
  % Dimensions
  f.Dimensions(1).Name      = 'NX';
  f.Dimensions(1).Length    = grid.nx;
  f.Dimensions(1).Unlimited = false;
  
  f.Dimensions(2).Name      = 'NY';
  f.Dimensions(2).Length    = grid.ny;
  f.Dimensions(2).Unlimited = false;
  
  f.Dimensions(3).Name      = 'time';
  f.Dimensions(3).Length    = 2;
  f.Dimensions(3).Unlimited = true;
  
  % Dimension variables
  
  % x
  f.Variables(1).Name         = 'x';
  f.Variables(1).Dimensions   = f.Dimensions(1);
  f.Variables(1).Size         = grid.nx;
  f.Variables(1).Datatype     = 'double';
  f.Variables(1).Attributes(1).Name  = 'long_name';
  f.Variables(1).Attributes(1).Value = 'X-coordinate';
  f.Variables(1).Attributes(2).Name  = 'units';
  f.Variables(1).Attributes(2).Value = 'm';
  f.Variables(1).ChunkSize    = [];
  f.Variables(1).FillValue    = [];
  f.Variables(1).DeflateLevel = [];
  f.Variables(1).Shuffle      = false;
  
  % y
  f.Variables(2).Name         = 'y';
  f.Variables(2).Dimensions   = f.Dimensions(2);
  f.Variables(2).Size         = grid.ny;
  f.Variables(2).Datatype     = 'double';
  f.Variables(2).Attributes(1).Name  = 'long_name';
  f.Variables(2).Attributes(1).Value = 'Y-coordinate';
  f.Variables(2).Attributes(2).Name  = 'units';
  f.Variables(2).Attributes(2).Value = 'm';
  f.Variables(2).ChunkSize    = [];
  f.Variables(2).FillValue    = [];
  f.Variables(2).DeflateLevel = [];
  f.Variables(2).Shuffle      = false;
  
  % time
  f.Variables(3).Name         = 'time';
  f.Variables(3).Dimensions   = f.Dimensions(3);
  f.Variables(3).Size         = 2;
  f.Variables(3).Datatype     = 'double';
  f.Variables(3).Attributes(1).Name  = 'long_name';
  f.Variables(3).Attributes(1).Value = 'Time';
  f.Variables(3).Attributes(2).Name  = 'units';
  f.Variables(3).Attributes(2).Value = 'years';
  f.Variables(3).ChunkSize    = [];
  f.Variables(3).FillValue    = [];
  f.Variables(3).DeflateLevel = [];
  f.Variables(3).Shuffle      = false;
  
  % SMB variable
  
  f.Variables(4).Name         = 'SMB';
  f.Variables(4).Dimensions   = f.Dimensions;
  f.Variables(4).Size         = [grid.nx, grid.ny, 2];
  f.Variables(4).Datatype     = 'double';
  f.Variables(4).Attributes(1).Name  = 'long_name';
  f.Variables(4).Attributes(1).Value = 'Surface mass balance';
  f.Variables(4).Attributes(2).Name  = 'units';
  f.Variables(4).Attributes(2).Value = 'mieq/yr';
  f.Variables(4).ChunkSize    = [];
  f.Variables(4).FillValue    = [];
  f.Variables(4).DeflateLevel = [];
  f.Variables(4).Shuffle      = false;
  
  % T2m variable
  
  f.Variables(5).Name         = 'T2m';
  f.Variables(5).Dimensions   = f.Dimensions;
  f.Variables(5).Size         = [grid.nx, grid.ny, 2];
  f.Variables(5).Datatype     = 'double';
  f.Variables(5).Attributes(1).Name  = 'long_name';
  f.Variables(5).Attributes(1).Value = '2-m annual mean temperature';
  f.Variables(5).Attributes(2).Name  = 'units';
  f.Variables(5).Attributes(2).Value = 'K';
  f.Variables(5).ChunkSize    = [];
  f.Variables(5).FillValue    = [];
  f.Variables(5).DeflateLevel = [];
  f.Variables(5).Shuffle      = false;
  
  % Final metadata
  f.Attributes = [];
  f.Groups     = [];
  f.Format     = 'classic';
  
end