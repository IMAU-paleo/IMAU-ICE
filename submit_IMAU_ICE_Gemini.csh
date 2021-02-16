rm -f run_IMAU_ICE_Gemini.sh.*

module load mpi/openmpi-x86_64

#qsub -cwd -m e ./run_IMAU_ICE_Gemini.sh
qsub -cwd -V ./run_IMAU_ICE_Gemini.sh