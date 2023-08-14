#! /bin/csh -f

# Create module and object folders if not already there (-p).
# Bonus: it creates also the src folder if it does not exist,
# although in that case you would have other, bigger problems.

mkdir -p src/module-files
mkdir -p src/object-files

# Go where the source code is

cd src

# Remove any existing compiled files

make clean

# Compile the code using settings in src/Makefile. The compiled
# files will be stored in the module and object folders.

make all

# Go back to main model directory

cd ..

# Remove any existing, older model executables

rm -f IMAU_ICE_program

# Move the newly compiled executable from the src folder
# into the main model directory

mv src/IMAU_ICE_program .

# Now you are ready to run the model using something like
# mpiexec -n <number_of_cores> IMAU_ICE_program <path_to_config_file>
# Have fun :)
