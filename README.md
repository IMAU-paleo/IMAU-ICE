# IMAU-ICE

Last updated: 2021-02-16 by Tijn Berends (c.j.berends@uu.nl)

This repository contains the source code of IMAU-ICE, as well as a few scripts for compiling and running the code locally and on the UU Gemini system.

===========================
===== GETTING STARTED =====
===========================

1.  Create folder called "IMAU-ICE" somewhere on your computer.
2.  Download the "src_v1.1.0"   folder from the GitHub repository and place it in the IMAU-ICE folder (so you have IMAU-ICE/src_v1.1.0/XXX.f90, etc.)
3.  Download the "config-files" folder from the GitHub repository and place it in the IMAU-ICE folder (so you have IMAU-ICE/config-files/config_test_small, etc.)
4.  Download the scripts "compile_clean.csh" and "run_IMAU_ICE_local.csh" and place the in the IMAU-ICE folder.
5.  Open IMAU-ICE/src_v1.1.0/Makefile using your favourite text editor.
6.  Make sure that Makefile.include-local is included (and all others are commented out).
7.  Run "compile_clean.csh" in the IMAU-ICE folder by typing "./compile_clean.csh" in the terminal.
8a. If it works: congratulations, you've succesfully compiled the model.
8b. If it doesn't: this most likely is because the external packages needed by the compiler (Lapack, NetCDF, and MPI) are located somewhere else in your system than they are in mine. In that case, ask your local tech-savvy IT person to fix this.
9.  Once the model has succesfully compiled, run it by typing "./run_IMAU_ICE_local.csh". If succesfull, you should see the following output: (see below).
10. An output folder should now have appeared in your IMAU-ICE folder (IMAU-ICE/results_YYYYMMDD_001) containing several ASCII text files and two NetCDF files.

This runs a short 1000-year schematic experiment. If all of this works, you can move on to realistic simulations. For this, you need appropriate input files describing present-day topography/mathymetry, ice geometry, and climate, as well as some sort of forcing record. Contact me if you need these.


Screen output for the test run in step 9:
 
 ====================================================
 ===== Running IMAU-ICE v1.1.0_dev on   2 cores =====
 ====================================================
 
  Running benchmark experiment "EISMINT_1"
 
  Output directory: results_20210216_002/
 
  Initialising insolation data from /Users/berends/Documents/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc...
 
  Initialising model region ANT (Antarctica)...
   Initialising model grid at 50.00 km resolution: [  29 x   29] pixels
   Mapping PD      data to model grid...
   Mapping init    data to model grid...
   Initialising climate model...
   Initialising SMB model...
   Initialising BMB model...
   Initialising ice dynamics model...
   Initialising ELRA GIA model...
  Finished initialising model region ANT.
 
 Coupling model: t =  -120.000 kyr
 
  Running model region ANT (   Antarctica) from t =  -120.000 to t =  -119.900 kyr
   t =  -120.00 kyr - writing output...
 
 Coupling model: t =  -119.900 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.900 to t =  -119.800 kyr
   t =  -119.90 kyr - writing output...
 
 Coupling model: t =  -119.800 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.800 to t =  -119.700 kyr
   t =  -119.80 kyr - writing output...
 
 Coupling model: t =  -119.700 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.700 to t =  -119.600 kyr
   t =  -119.70 kyr - writing output...
 
 Coupling model: t =  -119.600 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.600 to t =  -119.500 kyr
   t =  -119.60 kyr - writing output...
 
 Coupling model: t =  -119.500 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.500 to t =  -119.400 kyr
   t =  -119.50 kyr - writing output...
 
 Coupling model: t =  -119.400 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.400 to t =  -119.300 kyr
   t =  -119.40 kyr - writing output...
 
 Coupling model: t =  -119.300 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.300 to t =  -119.200 kyr
   t =  -119.30 kyr - writing output...
 
 Coupling model: t =  -119.200 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.200 to t =  -119.100 kyr
   t =  -119.20 kyr - writing output...
 
 Coupling model: t =  -119.100 kyr
 
  Running model region ANT (   Antarctica) from t =  -119.100 to t =  -119.000 kyr
   t =  -119.10 kyr - writing output...
   t =  -119.00 kyr - writing output...

 ================================================================================
 ===== Simulation finished in  0 days,  0 hours,  0 minutes and  2 seconds! =====
 ================================================================================
