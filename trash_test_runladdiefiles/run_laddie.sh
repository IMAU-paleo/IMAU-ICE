#!/bin/bash

# Activate conda environment
source /Users/5941962/opt/anaconda3/etc/profile.d/conda.sh
conda activate

# Check if help_fields has time dimension, if so: copy that one to imauice_output.nc'
echo '... Copying updated help_fields_ANT.nc to imauice_output.nc'
cd 
cd /Users/5941962/surfdrive/AA_COUPLING_V2
python3 check_for_helpfields.py '/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4/help_fields_ANT.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/imauice_output.nc'
wait

# Prepair current geometry for LADDIE:
echo '... Converting geometry to laddie_geom.nc ...'
cd 
cd /Users/5941962/surfdrive/AA_COUPLING_V2
python3 convert_geometry.py '/Users/5941962/surfdrive/AA_COUPLING_V2/imauice_output.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/laddie_geom.nc'
wait
echo '... Converting geometry to laddie_geom.nc finished!'


# Run LADDIE
echo 'Computing melt rates using LADDIE'
cd
cd /Users/5941962/surfdrive/laddie
echo 'Waiting for LADDIE ...'
python3 runladdie.py config_testv2.toml
wait
echo '... LADDIE is finished!'

input_dir="/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4/testv2/"
output_dir="/Users/5941962/surfdrive/AA_COUPLING_V2/"

# Find the file with the highest number
latest_file=$(ls "$input_dir"/output_*.nc | sort -V | tail -n 1)

# Check if any files were found
if [ -n "$latest_file" ]; then
  # Copy the latest file to the output directory
  cp "$latest_file" "$output_dir"/melt_raw.nc
  echo "Latest file '$latest_file' copied successfully."
else
  echo "No files found in '$input_dir'."
fi

# Find the file with the highest number
latest_file=$(ls "$input_dir"/restart_*.nc | sort -V | tail -n 1)

# Check if any files were found
if [ -n "$latest_file" ]; then
  # Copy the latest file to the output directory
  cp "$latest_file" "$output_dir"/laddie_restart.nc
  echo "Latest file '$latest_file' copied successfully."
else
  echo "No files found in '$input_dir'."
fi


# Prepair LADDIE output for IMAU-ICE (create laddie.nc)
echo '... Converting meltfield to laddie_melt.nc ...'
#cp /Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4/testv2/output_*.nc /Users/5941962/surfdrive/AA_COUPLING_V2/melt_raw.nc
cd 
cd /Users/5941962/surfdrive/AA_COUPLING_V2
python3 convert_meltfield.py '/Users/5941962/surfdrive/AA_COUPLING_V2/melt_raw.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/laddie_melt.nc'
wait
echo '... Converting meltfield to laddie_melt.nc finished!'

#cp /Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4/testv2/restart_*.nc /Users/5941962/surfdrive/AA_COUPLING_V2/laddie_restart.nc

#rm -r /Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv4/testv2

conda deactivate

# End of script