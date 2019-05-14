########  SCRIPT EXAMPLE ##################

# @ shell=/bin/bash
# @ error   = job100T1p5.err.$(jobid)
# @ output  = job100T1p5.out.$(jobid)
# @ job_type = parallel
# @ environment = COPY_ALL
# @ node_usage= shared
# @ node = 1
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit =1:00:00
# @ notification = complete
# @ node_resources = ConsumableMemory(56gb)
# @ queue

module load intel/14.0
module load python27/python/2.7
module load python27/scipy/2015.10
  
# set all ESPResSo++ environment variables

source /ptmp/gvargas/work/collabJSmrek/e++ShWalls2Therms2k16JuneInst/ESPRC

mpiexec -n 16 python 2tmf_v3Inst.py 1000 100 1.5 0.005 1 3320000

