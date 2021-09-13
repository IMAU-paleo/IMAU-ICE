#!/bin/csh

# Don't run this script by itself! Instead run "submit_IMAU_ICE_Cartesius.sh",
# that one calls this one with sbatch and the required input arguments,
# and it also loads the required modules.
                                                                   
#SBATCH -N 1
#SBATCH -t 00:59:00
set number_of_cores = 24

if(! $#argv == 2 ) then
 echo ' ERROR: needs two arguments, e.g.: sbatch run_IMAU_ICE_Cartesius.csh config_test name_of_run'
 exit 0
endif

# Check if the executable "IMAU_ICE_program" actually exists
if(! -e IMAU_ICE_program) then
 echo ' ERROR: executable "IMAU_ICE_program" doesnt exist!'
 exit 0
endif

# Config file and run name specified through input arguments
set config_file_name = $1
set run_name         = $2

# Check to make sure a results directory of this name doesn't already exist
set run_dir_name = 'results_IMAU_ICE_Cartesius_'${run_name}
if( -d $TMPDIR/${run_dir_name}) then
 echo ' WARNING: directory '${run_dir_name}' already exists on the scratch disk "' $TMPDIR '; deleting it!'
 rm -rf $TMPDIR/${run_dir_name}
 wait
endif
if( -d $HOME/IMAU_ICE/${run_dir_name}) then
 echo ' WARNING: directory '${run_dir_name}' already exists in your IMAU-ICE directory; deleting it!'
 rm -rf $HOME/IMAU_ICE/${run_dir_name}
 wait
endif

# Check if this config file actually exists
if(! -e ${config_file_name}) then
 echo ' ERROR: config file "'${config_file_name}'" doesnt exist!'
 exit 0
endif

# All required stuff seems to exist. Continue.

# Write some key info to the slurm-XXXXXX.out (output streams are automatically redirected here by SLURM, nothing to be done about that)
echo '==========================='
echo '===== IMAU_ICE on Cartesius ====='
echo '==========================='
echo ''
echo '  run name       : '${run_name}
echo '  config file    : '${config_file_name}
echo ''

echo '  Started on     : '
date +"%m-%d-%y, %T"
echo ''
echo ''

# Create a directory for this run on the scratch disk and copy the executable and config file there
echo '  Creating run directory on scratch disc '$TMPDIR'/'${run_dir_name}
mkdir $TMPDIR/${run_dir_name}
cp $HOME/IMAU_ICE/IMAU_ICE_program     $TMPDIR/${run_dir_name}/
cp $HOME/IMAU_ICE/${config_file_name}  $TMPDIR/${run_dir_name}/
chmod 777 $TMPDIR/${run_dir_name}/IMAU_ICE_program

# Create a directory for this run on the home disk and copy the executable and config file there
mkdir $HOME/IMAU_ICE/${run_dir_name}
cp $HOME/IMAU_ICE/IMAU_ICE_program     $HOME/IMAU_ICE/${run_dir_name}/
cp $HOME/IMAU_ICE/${config_file_name}  $HOME/IMAU_ICE/${run_dir_name}/

# Move to the run directory on the scratch disk
cd $TMPDIR/${run_dir_name}

# Run the program
srun -n ${number_of_cores} IMAU_ICE_program $HOME/IMAU_ICE/${config_file_name} 

# Copy output files back to the home results directory
cp results/*.txt   $HOME/IMAU_ICE/${run_dir_name}
cp results/*.nc    $HOME/IMAU_ICE/${run_dir_name}

# Remove the temporary run directory from the sratch disk
rm -rf $TMPDIR/${run_dir_name}

echo ''
echo '  Finished on    : '
date +"%m-%d-%y, %T"
echo ''

# Wait until everything is finished
wait
