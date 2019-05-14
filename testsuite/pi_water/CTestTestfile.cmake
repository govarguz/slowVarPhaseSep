# CMake generated Testfile for 
# Source directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pi_water
# Build directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pi_water
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pi_water "/usr/bin/python2" "/data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/pi_water/water.py")
set_tests_properties(pi_water PROPERTIES  ENVIRONMENT "PYTHONPATH=/data/isilon/vargas/phaseSep/e++ShWallsInst201606:/data/isilon/vargas/phaseSep/e++ShWallsInst201606/contrib:/home/theorie/vargas/espressopp2/espAsItWas2k15:/people/thnfs/homes/vargas/.local/lib/python2.7/site-packages")
