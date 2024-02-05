#!/bin/bash

#==============================================================================================================
# Activate conda environment !! Depends on user !! 
source /Users/5941962/opt/anaconda3/etc/profile.d/conda.sh
conda activate


#==============================================================================================================
# Specify run_names:
laddie_run_name=testv2

# Specify directories
COUPLING_folder=/Users/5941962/surfdrive/AA_COUPLING_V2
OUTPUT_folder_IMAUICE=/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4
OUTPUT_folder_laddie=/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4'/'$laddie_run_name
DIR_laddie=/Users/5941962/surfdrive/laddie

# Specify config file laddie
laddie_CONFIG=config_testv2.toml

#==============================================================================================================
# Check if help_fields has time dimension, if so: copy that one to imauice_output.nc'
echo '... Copying updated help_fields_ANT.nc to imauice_output.nc'
cd 
cd $COUPLING_folder #/Users/5941962/surfdrive/AA_COUPLING_V2
python3 check_for_helpfields.py $OUTPUT_folder_IMAUICE'/help_fields_ANT.nc' $COUPLING_folder'/imauice_output.nc'
wait


#==============================================================================================================
# Prepair current geometry for LADDIE:
echo '... Converting geometry to laddie_geom.nc ...'
cd 
cd $COUPLING_folder
python3 convert_geometry.py $COUPLING_folder'/imauice_output.nc' $COUPLING_folder'/laddie_geom.nc'
wait
echo '... Converting geometry to laddie_geom.nc finished!'


#==============================================================================================================
# Run LADDIE
echo 'Computing melt rates using LADDIE'
cd
cd $DIR_laddie
echo 'Waiting for LADDIE ...'
python3 runladdie.py $laddie_CONFIG
wait
echo '... LADDIE is finished!'


#==============================================================================================================
# Copy output and restart file to coupling folder
# Find the file with the highest number
latest_file=$(ls "$OUTPUT_folder_laddie"/output_*.nc | sort -V | tail -n 1)

# Check if any files were found
if [ -n "$latest_file" ]; then
  # Copy the latest file to the output directory
  cp "$latest_file" "$COUPLING_folder"/melt_raw.nc
  echo "Latest file '$latest_file' copied successfully."
else
  echo "No files found in '$OUTPUT_folder_laddie'."
fi

# Find the file with the highest number
latest_file=$(ls "$OUTPUT_folder_laddie"/restart_*.nc | sort -V | tail -n 1)

# Check if any files were found
if [ -n "$latest_file" ]; then
  # Copy the latest file to the output directory
  cp "$latest_file" "$COUPLING_folder"/laddie_restart.nc
  echo "Latest file '$latest_file' copied successfully."
else
  echo "No files found in '$OUTPUT_folder_laddie'."
fi


#==============================================================================================================
# Prepair LADDIE output for IMAU-ICE (create laddie.nc)
echo '... Converting meltfield to laddie_melt.nc ...'
cd 
cd $COUPLING_folder
python3 convert_meltfield.py $COUPLING_folder'/melt_raw.nc' $COUPLING_folder'/laddie_melt.nc'
wait
echo '... Converting meltfield to laddie_melt.nc finished!'


conda deactivate

# End of script