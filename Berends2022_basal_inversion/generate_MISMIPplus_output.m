function generate_MISMIPplus_output( foldername)

% clc
% clear all
% close all
% 
% foldername = 'MISMIPplus_ice0_2km_slid1';

disp(['  Processing MISMIP+ output for run ' foldername '...'])

%% Read data
filename = [foldername '/restart_ANT.nc'];

x    = ncread(filename,'x'); x = x - min(x); dx = abs(x(2)-x(1));
y    = ncread(filename,'y');
time = ncread(filename,'time');

iceVolume    = zeros( length(time),1);
iceVAF       = zeros( length(time),1);
groundedArea = zeros( length(time),1);

xGL = zeros( length(time),length(y));
yGL = zeros( length(time),length(y));

for ti = 1:length(time)
  
  % Read data
  Hi = ncread(filename,'Hi',[1,1,ti],[Inf,Inf,1]);
  Hb = ncread(filename,'Hb',[1,1,ti],[Inf,Inf,1]);
  TAF = Hi - max(0, (0 - Hb) * (1028 / 910));
  
  % Scalars
  iceVolume(    ti) = sum(Hi(:)) * dx^2;
  iceVAF(       ti) = sum(TAF(TAF>0)) * dx^2;
  groundedArea( ti) = sum(TAF(:)>0) * dx^2;
  
  % GL
  yGL( ti,:) = y;
  for j = 1:length(y)
    i2 = 1;
    while TAF(i2,j)>0; i2 = i2+1; end
    i1 = i2-1;
    if i1==0
      disp('beep')
    end
    lambda = TAF(i1,j) / (TAF(i1,j)-TAF(i2,j));
    xGL( ti,j) = lambda * x(i2) + (1-lambda) * x(i1);
  end
  
end

%% Create NetCDF template
f = ncinfo('/Users/berends/Documents/Datasets/MISMIP+/supplement/submission_data/Ice0_CBO_ISSM_SSA_Tsai_500m.nc');

% nPointGL dimension
f.Dimensions(1).Length = length(y);
% Time dimension
f.Dimensions(2).Length = length(time);

% Time variable
f.Variables(1).Dimensions = f.Dimensions(2);
f.Variables(1).Size = length(time);

% Scalar variables (iceVolume, iceVAF, groundedArea)
f.Variables(2).Dimensions = f.Dimensions(2);
f.Variables(2).Size = length(time);
f.Variables(3).Dimensions = f.Dimensions(2);
f.Variables(3).Size = length(time);
f.Variables(4).Dimensions = f.Dimensions(2);
f.Variables(4).Size = length(time);

% xGL
f.Variables(5).Dimensions(1) = f.Dimensions(2);
f.Variables(5).Dimensions(2) = f.Dimensions(1);
f.Variables(5).Size = [length(time), length(y)];

% yGL
f.Variables(6).Dimensions(1) = f.Dimensions(2);
f.Variables(6).Dimensions(2) = f.Dimensions(1);
f.Variables(6).Size = [length(time), length(y)];

%% Create file and write data
filename = [foldername '/MISMIPplus_output.nc'];

if exist(filename,'file')
  delete(filename)
end

ncwriteschema(filename,f);

ncwrite(filename,'time'        ,time        );
ncwrite(filename,'iceVolume'   ,iceVolume   );
ncwrite(filename,'iceVAF'      ,iceVAF      );
ncwrite(filename,'groundedArea',groundedArea);
ncwrite(filename,'xGL'         ,xGL         );
ncwrite(filename,'yGL'         ,yGL         );

end