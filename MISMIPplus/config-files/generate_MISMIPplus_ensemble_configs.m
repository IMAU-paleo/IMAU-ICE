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

ensemble_name = 'MISMIPplus';

%% Different resolution/sliding law/stress balance/flow factor combinations for MISMIP+
Dims{ 1}.Opts{ 1}.suffix         = '_5km_slid1_DIVA';
Dims{ 1}.Opts{ 1}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 1}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 1}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 1}.Vars{ 2}.Value = '1';
Dims{ 1}.Opts{ 1}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 1}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 1}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 1}.Vars{ 4}.Value = '1.262E-17';
Dims{ 1}.Opts{ 1}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 1}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid1_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 2}.suffix         = '_5km_slid2_DIVA';
Dims{ 1}.Opts{ 2}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 2}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 2}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 2}.Vars{ 2}.Value = '2';
Dims{ 1}.Opts{ 2}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 2}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 2}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 2}.Vars{ 4}.Value = '1.038E-17';
Dims{ 1}.Opts{ 2}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 2}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid2_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 3}.suffix         = '_5km_slid3_DIVA';
Dims{ 1}.Opts{ 3}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 3}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 3}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 3}.Vars{ 2}.Value = '3';
Dims{ 1}.Opts{ 3}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 3}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 3}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 3}.Vars{ 4}.Value = '1.034E-17';
Dims{ 1}.Opts{ 3}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 3}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid3_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 4}.suffix         = '_2km_slid1_DIVA';
Dims{ 1}.Opts{ 4}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 4}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 4}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 4}.Vars{ 2}.Value = '1';
Dims{ 1}.Opts{ 4}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 4}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 4}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 4}.Vars{ 4}.Value = '1.517E-17';
Dims{ 1}.Opts{ 4}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 4}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid1_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 5}.suffix         = '_2km_slid2_DIVA';
Dims{ 1}.Opts{ 5}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 5}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 5}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 5}.Vars{ 2}.Value = '2';
Dims{ 1}.Opts{ 5}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 5}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 5}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 5}.Vars{ 4}.Value = '1.366E-17';
Dims{ 1}.Opts{ 5}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 5}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid2_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 6}.suffix         = '_2km_slid3_DIVA';
Dims{ 1}.Opts{ 6}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 6}.Vars{ 1}.Value = '2000';
Dims{ 1}.Opts{ 6}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 6}.Vars{ 2}.Value = '3';
Dims{ 1}.Opts{ 6}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 6}.Vars{ 3}.Value = 'DIVA';
Dims{ 1}.Opts{ 6}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 6}.Vars{ 4}.Value = '1.375E-17';
Dims{ 1}.Opts{ 6}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 6}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid3_DIVA/restart_ANT.nc';

Dims{ 1}.Opts{ 7}.suffix         = '_5km_slid1_SIASSA';
Dims{ 1}.Opts{ 7}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 7}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 7}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 7}.Vars{ 2}.Value = '1';
Dims{ 1}.Opts{ 7}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 7}.Vars{ 3}.Value = 'SIA/SSA';
Dims{ 1}.Opts{ 7}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 7}.Vars{ 4}.Value = '1.211E-17';
Dims{ 1}.Opts{ 7}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 7}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid1_SIASSA/restart_ANT.nc';

Dims{ 1}.Opts{ 8}.suffix         = '_5km_slid2_SIASSA';
Dims{ 1}.Opts{ 8}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 8}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 8}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 8}.Vars{ 2}.Value = '2';
Dims{ 1}.Opts{ 8}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 8}.Vars{ 3}.Value = 'SIA/SSA';
Dims{ 1}.Opts{ 8}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 8}.Vars{ 4}.Value = '9.962E-18';
Dims{ 1}.Opts{ 8}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 8}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid2_SIASSA/restart_ANT.nc';

Dims{ 1}.Opts{ 9}.suffix         = '_5km_slid3_SIASSA';
Dims{ 1}.Opts{ 9}.Vars{ 1}.Name  = 'dx_ANT_config';
Dims{ 1}.Opts{ 9}.Vars{ 1}.Value = '5000';
Dims{ 1}.Opts{ 9}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
Dims{ 1}.Opts{ 9}.Vars{ 2}.Value = '3';
Dims{ 1}.Opts{ 9}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
Dims{ 1}.Opts{ 9}.Vars{ 3}.Value = 'SIA/SSA';
Dims{ 1}.Opts{ 9}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
Dims{ 1}.Opts{ 9}.Vars{ 4}.Value = '9.898E-18';
Dims{ 1}.Opts{ 9}.Vars{ 5}.Name  = 'filename_init_ANT_config';
Dims{ 1}.Opts{ 9}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_5km_slid3_SIASSA/restart_ANT.nc';

% Dims{ 1}.Opts{10}.suffix         = '_2km_slid1_SIASSA';
% Dims{ 1}.Opts{10}.Vars{ 1}.Name  = 'dx_ANT_config';
% Dims{ 1}.Opts{10}.Vars{ 1}.Value = '2000';
% Dims{ 1}.Opts{10}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
% Dims{ 1}.Opts{10}.Vars{ 2}.Value = '1';
% Dims{ 1}.Opts{10}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
% Dims{ 1}.Opts{10}.Vars{ 3}.Value = 'SIA/SSA';
% Dims{ 1}.Opts{10}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
% Dims{ 1}.Opts{10}.Vars{ 4}.Value = '1.517E-17';
% Dims{ 1}.Opts{10}.Vars{ 5}.Name  = 'filename_init_ANT_config';
% Dims{ 1}.Opts{10}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid1_SIASSA/restart_ANT.nc';
% 
% Dims{ 1}.Opts{11}.suffix         = '_2km_slid2_SIASSA';
% Dims{ 1}.Opts{11}.Vars{ 1}.Name  = 'dx_ANT_config';
% Dims{ 1}.Opts{11}.Vars{ 1}.Value = '2000';
% Dims{ 1}.Opts{11}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
% Dims{ 1}.Opts{11}.Vars{ 2}.Value = '2';
% Dims{ 1}.Opts{11}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
% Dims{ 1}.Opts{11}.Vars{ 3}.Value = 'SIA/SSA';
% Dims{ 1}.Opts{11}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
% Dims{ 1}.Opts{11}.Vars{ 4}.Value = '1.366E-17';
% Dims{ 1}.Opts{11}.Vars{ 5}.Name  = 'filename_init_ANT_config';
% Dims{ 1}.Opts{11}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid2_SIASSA/restart_ANT.nc';
% 
% Dims{ 1}.Opts{12}.suffix         = '_2km_slid3_SIASSA';
% Dims{ 1}.Opts{12}.Vars{ 1}.Name  = 'dx_ANT_config';
% Dims{ 1}.Opts{12}.Vars{ 1}.Value = '2000';
% Dims{ 1}.Opts{12}.Vars{ 2}.Name  = 'MISMIPplus_sliding_law_config';
% Dims{ 1}.Opts{12}.Vars{ 2}.Value = '3';
% Dims{ 1}.Opts{12}.Vars{ 3}.Name  = 'choice_ice_dynamics_config';
% Dims{ 1}.Opts{12}.Vars{ 3}.Value = 'SIA/SSA';
% Dims{ 1}.Opts{12}.Vars{ 4}.Name  = 'MISMIPplus_A_flow_initial_config';
% Dims{ 1}.Opts{12}.Vars{ 4}.Value = '1.375E-17';
% Dims{ 1}.Opts{12}.Vars{ 5}.Name  = 'filename_init_ANT_config';
% Dims{ 1}.Opts{12}.Vars{ 5}.Value = 'MISMIPplus/MISMIPplus_init_2km_slid3_SIASSA/restart_ANT.nc';

%% Different sub-grid melt schemes
Dims{ 2}.Opts{ 1}.suffix         = '_NMP';
Dims{ 2}.Opts{ 1}.Vars{ 1}.Name  = 'choice_BMB_subgrid_config';
Dims{ 2}.Opts{ 1}.Vars{ 1}.Value = 'NMP';

Dims{ 2}.Opts{ 2}.suffix         = '_FCMP';
Dims{ 2}.Opts{ 2}.Vars{ 1}.Name  = 'choice_BMB_subgrid_config';
Dims{ 2}.Opts{ 2}.Vars{ 1}.Value = 'FCMP';

Dims{ 2}.Opts{ 3}.suffix         = '_PMP';
Dims{ 2}.Opts{ 3}.Vars{ 1}.Name  = 'choice_BMB_subgrid_config';
Dims{ 2}.Opts{ 3}.Vars{ 1}.Value = 'PMP';

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