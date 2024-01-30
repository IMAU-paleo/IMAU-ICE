#!/bin/bash

# Activate conda environment
source /Users/5941962/opt/anaconda3/etc/profile.d/conda.sh
conda activate

# Copy help_fields_ANT.nc 
if [ -e  '/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv3/help_fields_ANT.nc' ]
	then cp '/Users/5941962/surfdrive/IMAU-ICE/test_LADDIEv3/help_fields_ANT.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/imauice_output.nc'
fi

# Prepair current geometry for LADDIE:
echo '... Converting geometry to laddie_geom.nc ...'
cd 
cd /Users/5941962/surfdrive/AA_COUPLING_V2
python3 convert_geometry.py '/Users/5941962/surfdrive/AA_COUPLING_V2/imauice_output.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/laddie_geom.nc'
wait
echo '... Converting geometry to laddie_geom.nc finished!'


# Prepair config file LADDIE (depends on forcing)


# Run LADDIE
echo 'Computing melt rates using LADDIE'
cd
cd /Users/5941962/surfdrive/laddie
echo 'Waiting for LADDIE ...'
python3 runladdie.py config_testv2.toml
wait
echo '... LADDIE is finished!'



# Prepair LADDIE output for IMAU-ICE (create laddie.nc)
echo '... Converting meltfield to laddie_melt.nc ...'
cp /Users/5941962/surfdrive/laddie/output/testv2/output_*.nc /Users/5941962/surfdrive/AA_COUPLING_V2/melt_raw.nc
cd 
cd /Users/5941962/surfdrive/AA_COUPLING_V2
python3 convert_meltfield.py '/Users/5941962/surfdrive/AA_COUPLING_V2/melt_raw.nc' '/Users/5941962/surfdrive/AA_COUPLING_V2/laddie_melt.nc'
wait
echo '... Converting meltfield to laddie_melt.nc finished!'

cp /Users/5941962/surfdrive/laddie/output/testv2/restart_*.nc /Users/5941962/surfdrive/AA_COUPLING_V2/laddie_restart.nc
# 

rm -r /Users/5941962/surfdrive/laddie/output/testv2

conda deactivate

# End of script