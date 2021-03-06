#!/bin/sh
#$ -cwd
#$ -pe PE_Multi_IB 320
#$ -N isPolGRZG1
#$ -j y
#$ -o logH.$JOB_ID.out
#$ -m e
# #$ -M vargas@mpip-mainz.mpg.de
#$ -l h_rt=36:00:00
#$ -S /bin/bash
source /data/isilon/vargas/adressIB1/ESPRC
source /sw/linux/modules/init/bash
module load openmpi
/sw/linux/mpi/gcc/openmpi/bin/mpirun -np 320 -x PYTHONPATH -x PATH -x LD_LIBRARY_PATH python polymer_meltPoC320.py
