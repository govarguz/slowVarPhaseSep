########  SCRIPT EXAMPLE ##################

# @ shell=/bin/bash
# @ error   = job.err.$(jobid)
# @ output  = job.out.$(jobid)
# @ job_type = parallel
# @ environment = COPY_ALL
# @ node_usage= shared
# @ node = 1
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit = 01:10:00
# @ notification = complete
# @ queue

module load intel/14.0
module load python27/python/2.7
  
# set all ESPResSo++ environment variables

source /ptmp/gvargas/work/collabJSmrek/e++ShWalls2Therms2k16June/ESPRC

mpiexec -n 16 python 2tmf_v2.py

