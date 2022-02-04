clc
clear all
close all

%% Set up the ensemble
%
% An ensemble has different "dimensions", each of which has several
% "options". The total number of ensemble members is therefore equal to
% the sum of all options over all dimensions.
% An option can contain several variables, because often you'll want to
% change several variables at the same time (e.g. resolution, sliding law, and flow
% factor in the MISMIP+ experiments)

ensemble_name = 'MISOMIP1';

%% Different resolution/sliding law/flow factor combinations for MISMIP+
Dims{ 1}.Opts{ 1}.suffix         = '_5km_slid1';
Dims{ 1}.Opts{ 1}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 1}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 1}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 1}.Vars{ 2}.Value = '1';
Dims{ 1}.Opts{ 1}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 1}.Vars{ 3}.Value = '1.262E-17';
Dims{ 1}.Opts{ 1}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 1}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid1/restart_ANT.nc';

Dims{ 1}.Opts{ 2}.suffix         = '_5km_slid2';
Dims{ 1}.Opts{ 2}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 2}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 2}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 2}.Vars{ 2}.Value = '2';
Dims{ 1}.Opts{ 2}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 2}.Vars{ 3}.Value = '1.038E-17';
Dims{ 1}.Opts{ 2}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 2}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid2/restart_ANT.nc';

Dims{ 1}.Opts{ 3}.suffix         = '_5km_slid3';
Dims{ 1}.Opts{ 3}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 3}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 3}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 3}.Vars{ 2}.Value = '3';
Dims{ 1}.Opts{ 3}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 3}.Vars{ 3}.Value = '1.034E-17';
Dims{ 1}.Opts{ 3}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 3}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid3/restart_ANT.nc';

Dims{ 1}.Opts{ 4}.suffix         = '_2km_slid1';
Dims{ 1}.Opts{ 4}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 4}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 4}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 4}.Vars{ 2}.Value = '1';
Dims{ 1}.Opts{ 4}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 4}.Vars{ 3}.Value = '1.517E-17';
Dims{ 1}.Opts{ 4}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 4}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid1/restart_ANT.nc';

Dims{ 1}.Opts{ 5}.suffix         = '_2km_slid2';
Dims{ 1}.Opts{ 5}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 5}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 5}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 5}.Vars{ 2}.Value = '2';
Dims{ 1}.Opts{ 5}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 5}.Vars{ 3}.Value = '1.366E-17';
Dims{ 1}.Opts{ 5}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 5}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid2/restart_ANT.nc';

Dims{ 1}.Opts{ 6}.suffix         = '_2km_slid3';
Dims{ 1}.Opts{ 6}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 6}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 6}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 6}.Vars{ 2}.Value = '3';
Dims{ 1}.Opts{ 6}.Vars{ 3}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 6}.Vars{ 3}.Value = '1.375E-17';
Dims{ 1}.Opts{ 6}.Vars{ 4}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 6}.Vars{ 4}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid3/restart_ANT.nc';

%% All the BMB parameterisations
Dims{ 2}.Opts{ 1}.suffix         = '_lin';
Dims{ 2}.Opts{ 1}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 1}.Vars{ 1}.Value = 'Favier2019_lin';
Dims{ 2}.Opts{ 1}.Vars{ 2}.Name  = 'BMB_Favier2019_lin_GammaT_config';
Dims{ 2}.Opts{ 1}.Vars{ 2}.Value = '3.3314E-05';

Dims{ 2}.Opts{ 2}.suffix         = '_quad';
Dims{ 2}.Opts{ 2}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 2}.Vars{ 1}.Value = 'Favier2019_quad';
Dims{ 2}.Opts{ 2}.Vars{ 2}.Name  = 'BMB_Favier2019_quad_GammaT_config';
Dims{ 2}.Opts{ 2}.Vars{ 2}.Value = '111.6E-5';

Dims{ 2}.Opts{ 3}.suffix         = '_Mplus';
Dims{ 2}.Opts{ 3}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 3}.Vars{ 1}.Value = 'Favier2019_Mplus';
Dims{ 2}.Opts{ 3}.Vars{ 2}.Name  = 'BMB_Favier2019_Mplus_GammaT_config';
Dims{ 2}.Opts{ 3}.Vars{ 2}.Value = '108.6E-5';

Dims{ 2}.Opts{ 4}.suffix         = '_plume';
Dims{ 2}.Opts{ 4}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 4}.Vars{ 1}.Value = 'Lazeroms2018_plume';
Dims{ 2}.Opts{ 4}.Vars{ 2}.Name  = 'BMB_Lazeroms2018_GammaT_config';
Dims{ 2}.Opts{ 4}.Vars{ 2}.Value = '5.493E-04';

Dims{ 2}.Opts{ 5}.suffix         = '_PICO';
Dims{ 2}.Opts{ 5}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 5}.Vars{ 1}.Value = 'PICO';
Dims{ 2}.Opts{ 5}.Vars{ 2}.Name  = 'BMB_PICO_GammaTstar_config';
Dims{ 2}.Opts{ 5}.Vars{ 2}.Value = '3.6131E-05';

Dims{ 2}.Opts{ 6}.suffix         = '_PICOP';
Dims{ 2}.Opts{ 6}.Vars{ 1}.Name  = 'choice_BMB_shelf_model_config';
Dims{ 2}.Opts{ 6}.Vars{ 1}.Value = 'PICOP';
Dims{ 2}.Opts{ 6}.Vars{ 2}.Name  = 'BMB_Lazeroms2018_GammaT_config';
Dims{ 2}.Opts{ 6}.Vars{ 2}.Value = '9.660E-04';

%% Clean up existing variation files
henk = dir;
namestart = [ensemble_name '_var'];
for i = 1:length(henk)
  if length(henk(i).name) > length(namestart)
    if strcmpi(henk(i).name(1:length(namestart)),namestart)
      delete(henk(i).name)
    end
  end
end

%% Generate the variation config files
generate_variation_files( Dims, ensemble_name)