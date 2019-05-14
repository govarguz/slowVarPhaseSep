#!/bin/sh
#$ -cwd
#$ -pe PE_Multi_IB 80
#$ -N ma3PolGRZG1
#$ -j y
#$ -o logH.$JOB_ID.out
#$ -m e
# #$ -M vargas@mpip-mainz.mpg.de
#$ -l h_rt=36:00:00
#$ -S /bin/bash
module purge
source /home/theorie/vargas/espressopp2/espMPItest62k15/ESPRC
source /sw/linux/modules/init/bash
module load openmpi
/sw/linux/mpi/gcc/openmpi/bin/mpirun -np 80 -x PYTHONPATH -x PATH -x LD_LIBRARY_PATH python polymer_meltPoC80.py
