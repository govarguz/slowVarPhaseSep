# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:32:40 2016

@author: vargas
"""
import numpy as np
import sys

cases = int(sys.argv[1])  # the folder name should start in 1
cores = int(sys.argv[2])  # In this case 216
nsteps = int(sys.argv[3])  # In this case is 10
neighbors=6 # Hardwire does not normaly change

###############################################################################
##############             Extracting NIntPerRank.dlb        ##################
###############################################################################
for k in range(1,cases+1):
	NInte1p=np.genfromtxt('T1p'+str(k)+'/NIntPerRank.dlb')
	Inte1pMax=[]
	[Inte1pMax.append(np.max(abs(NInte1p[i+1]-NInte1p[i]))) for i in range(nsteps-1)]
	print "InteMax1p",k,"es:", np.average(Inte1pMax)
	Inte1pAvg=[]
	[Inte1pAvg.append(np.average(abs(NInte1p[i+1]-NInte1p[i]))) for i in range(nsteps-1)]
	print "InteAvg1p",k,"es:", np.average(Inte1pAvg)

###############################################################################
##############             Extracting NPartsPerRank.dlb        ##################
###############################################################################

for k in range(1,cases+1):
	NParts1p=np.genfromtxt('T1p'+str(k)+'/NPartsPerRank.dlb')
	Parts1pMax=[]
	[Parts1pMax.append(np.max(abs(NParts1p[i+1]-NParts1p[i]))) for i in range(nsteps-1)]
	print "PartsMax1p",k,"es:", np.average(Parts1pMax)
	Parts1pAvg=[]
	[Parts1pAvg.append(np.average(abs(NParts1p[i+1]-NParts1p[i]))) for i in range(nsteps-1)]
	print "PartsAvg1p",k,"es:", np.average(Parts1pAvg)

###############################################################################
##############    Extracting PartsCommCellsPerRankPerDir.dlb ##################
###############################################################################

for k in range(1,cases+1):
	Comms1p= np.genfromtxt('T1p'+str(k)+'/PartsCommCellsPerRankPerDir.dlb',dtype=None,delimiter='"',deletechars="()",replace_space='',usecols=(range(1,cores*2,2)))
	Comms=np.empty([nsteps,cores,neighbors])
	for ii in range(nsteps):
	  for j in range(cores):
	    a=tuple(int(x) for x in Comms1p[ii][j][1:-1].split(','))
	    Comms[ii][j][:]=np.asarray(a)
	  print "done"
	Comms1pMax=[]
	[Comms1pMax.append(np.max(abs(Comms[i+1]-Comms[i]))) for i in range(9)]
	print "CommsMax1p",k,"es:", np.average(Comms1pMax)
	Comms1pAvg=[]
	[Comms1pAvg.append(np.average(abs(Comms[i+1]-Comms[i]))) for i in range(9)]
	print "CommsAvg1p",k,"es:", np.average(Comms1pAvg)

###############################################################################
##############    Extracting PartsHangoverPerRankPerDir.dlb ##################
###############################################################################

for k in range(1,cases+1):
	Hangs1p= np.genfromtxt('T1p'+str(k)+'/PartsHangoverPerRankPerDir.dlb',dtype=None,delimiter='"',deletechars="()",replace_space='',usecols=(range(1,216*2,2)))
	Hangs=np.empty([nsteps,cores,neighbors])
	for ii in range(nsteps):
	  for j in range(cores):
	    a=tuple(int(x) for x in Hangs1p[ii][j][1:-1].split(','))
	    Hangs[ii][j][:]=np.asarray(a)
	  print "done"
	Hangs1pMax=[]
	[Hangs1pMax.append(np.max(abs(Hangs[i+1]-Hangs[i]))) for i in range(1,9)]
	print "CommsMax1p",k,"es:", np.average(Hangs1pMax)
	Hangs1pAvg=[]
	[Hangs1pAvg.append(np.average(abs(Hangs[i+1]-Hangs[i]))) for i in range(1,9)]
	print "CommsAvg1p",k,"es:", np.average(Hangs1pAvg)

print "Done HAVG."