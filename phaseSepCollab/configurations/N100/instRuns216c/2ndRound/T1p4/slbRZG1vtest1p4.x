########  SCRIPT EXAMPLE ##################

# @ shell=/bin/bash
# @ error   = job100T1p4.err.$(jobid)
# @ output  = job100T1p4.out.$(jobid)
# @ job_type = parallel
# @ environment = COPY_ALL
# @ node_usage= not_shared
# @ node = 12
# @ tasks_per_node = 18
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,shared,us
# @ wall_clock_limit = 24:00:00
# @ notification = complete
# @ node_resources = ConsumableMemory(56gb)
# @ node_topology = island
# @ island_count = 1
# @ class = test
# @ queue

module load intel/14.0
module load python27/python/2.7
module load python27/scipy/2015.10
  
# set all ESPResSo++ environment variables

source /ptmp/gvargas/work/collabJSmrek/e++ShWalls2Therms2k16JuneInst/ESPRC

mpiexec -n 216 python 2tmf_v2Inst.py 1000 100 1.4 0.005 1 4410000

