#    Makefile for compiling the IMAU_ICE odel

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

#include Makefile.include-local
include Makefile.include-gemini

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ================================================================
#       Building the IMAU_ICE_program
# ================================================================

# Adding this program to the list of all programs:
all_programs += IMAU_ICE_program

# List of object files for fortran_example_model_using_namelist_program:
IMAU_ICE_program_objects = \
        configuration_module.o \
        parallel_module.o \
        parameters_module.o \
        data_types_netcdf_module.o \
        data_types_module.o \
        netcdf_module.o \
        forcing_module.o \
        global_text_output_module.o \
        utilities_module.o \
        derivatives_and_grids_module.o \
        reference_fields_module.o \
        general_ice_model_data_module.o \
        ice_dynamics_module.o \
        thermodynamics_module.o \
        SMB_module.o \
        BMB_module.o \
        climate_module.o \
        isotopes_module.o \
        bedrock_ELRA_module.o \
        IMAU_ICE_main_model.o \
        IMAU_ICE_program.o

# The executable (the target) depends on the object files:
IMAU_ICE_program : $(IMAU_ICE_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)
