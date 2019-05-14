# CMake generated Testfile for 
# Source directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pickle_potential
# Build directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pickle_potential
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pickle_potential "/usr/bin/python2" "/data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pickle_potential/pickle_potential.py")
set_tests_properties(pickle_potential PROPERTIES  ENVIRONMENT "PYTHONPATH=/data/isilon/vargas/phaseSep/e++ShWallsInst201606:/data/isilon/vargas/phaseSep/e++ShWallsInst201606/contrib:/home/theorie/vargas/espressopp2/espAsItWas2k15:/people/thnfs/homes/vargas/.local/lib/python2.7/site-packages")
add_test(testwarmup "/usr/bin/python2" "/data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pickle_potential/testwarmup.py")
set_tests_properties(testwarmup PROPERTIES  ENVIRONMENT "PYTHONPATH=/data/isilon/vargas/phaseSep/e++ShWallsInst201606:/data/isilon/vargas/phaseSep/e++ShWallsInst201606/contrib:/home/theorie/vargas/espressopp2/espAsItWas2k15:/people/thnfs/homes/vargas/.local/lib/python2.7/site-packages")
