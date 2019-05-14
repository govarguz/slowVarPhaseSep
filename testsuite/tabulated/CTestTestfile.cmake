# CMake generated Testfile for 
# Source directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/tabulated
# Build directory: /data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/tabulated
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(polymer_melt_tabulated "/usr/bin/python2" "/data/isilon/vargas/phaseSep/e++ShWallsInst201606/testsuite/tabulated/polymer_melt_tabulated.py")
set_tests_properties(polymer_melt_tabulated PROPERTIES  ENVIRONMENT "PYTHONPATH=/data/isilon/vargas/phaseSep/e++ShWallsInst201606:/data/isilon/vargas/phaseSep/e++ShWallsInst201606/contrib:/home/theorie/vargas/espressopp2/espAsItWas2k15:/people/thnfs/homes/vargas/.local/lib/python2.7/site-packages")
