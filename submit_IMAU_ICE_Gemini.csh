rm -f run_IMAU_ICE_Gemini.sh.*

module load mpi/openmpi-x86_64

#qsub -cwd -m e -V ./run_IMAU_ICE_Gemini.sh
#qsub -cwd -V ./run_IMAU_ICE_local.sh
qsub -cwd -m e -o ./run8.o -e run8.e -V ./run_IMAU_ICE_Gemini.sh
#qsub -cwd -m e -o ./run2.o -e run2.e -V -pe openmpi 16 ./run_IMAU_ICE_Gemini.sh 
