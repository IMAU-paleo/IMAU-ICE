# IMAU-ICE

Last updated: 2021-02-16 by Tijn Berends

This repository contains the source code of IMAU-ICE, as well as a few scripts for compiling and running the code locally and on the UU Gemini systems. In order for the code to succesfully compile, it requires access to the LAPACK, NetCDF, and MPI modules. This is arranged in the Makefile.include-XXX. On clusters like LISA and Gemini, you will need to load all of these modules plus your fortran compiler before compiling, and load MPI before running the model (simply typing ./program won't work with MPI programs).
