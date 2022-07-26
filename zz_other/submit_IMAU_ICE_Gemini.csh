#rm -f run_IMAU_ICE_Gemini.sh.*

module load mpi/openmpi-x86_64
module load petsc/3.16.3

#qsub -cwd -m e -V ./run_IMAU_ICE_Gemini.sh
#qsub -cwd -V ./run_IMAU_ICE_Gemini.sh

qsub -q all.q -cwd -m e -o ./run1.e -e ./run1.o -V ./run_IMAU_ICE_Gemini_1.sh
