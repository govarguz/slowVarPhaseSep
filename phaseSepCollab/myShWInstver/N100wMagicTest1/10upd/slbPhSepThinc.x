#!/bin/sh
#$ -cwd
#$ -pe PE_Multi_IB 64
#$ -N aPhSepb64E++2
#$ -j y
#$ -o logUJH.$JOB_ID.out
#$ -e logUJH.$JOB_ID.err
#$ -m e
# #$ -M vargas@mpip-mainz.mpg.de
#$ -l h_rt=36:00:00
#$ -S /bin/bash
source /data/isilon/vargas/phaseSep/e++ShWallsInst201606/ESPRC
source /sw/linux/modules/init/bash
module load openmpi
/sw/linux/mpi/gcc/openmpi/bin/mpirun -np 64 -x PYTHONPATH -x PATH -x LD_LIBRARY_PATH python 2tmf_v4Inst.py 1000 100 1.5 0.005 1 3320000 1

~                                                                                                        
