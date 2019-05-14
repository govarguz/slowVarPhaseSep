# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 13:21:19 2016

@author: vargas
"""
import espressopp
#import "_espressopp.so"
import mpi4py.MPI as MPI
import sys
import os
import time
from espressopp.tools.decomp import neiListHom
import csv

def myGatherConfOTF(system,velocities=True,unfolded=True,append=False,mpc=1):
  pidf=[]
  xposf=[]
  yposf=[]
  zposf=[]
  xvelf=[] 
  yvelf=[]
  zvelf=[]
  configurations = espressopp.analysis.ConfigurationsExt(system)
  configurations.unfolded = unfolded
  configurations.gather()
  configuration = configurations[0]

  if velocities:
    velocities = espressopp.analysis.Velocities(system)
    velocities.gather()
    velocity = velocities[0]

  numParticles  = int(espressopp.analysis.NPart(system).compute())
  box_x = system.bc.boxL[0]
  box_y = system.bc.boxL[1]
  box_z = system.bc.boxL[2]
  #st = "%d\n%15.10f %15.10f %15.10f\n" % (numParticles, box_x, box_y, box_z)
  #file.write(st)

  for pid in configuration:
        xpos   = configuration[pid][0]
        ypos   = configuration[pid][1]
        zpos   = configuration[pid][2]
        xposf.append(xpos)
	yposf.append(ypos)
	zposf.append(zpos)
        if velocities:
          xvel   = velocity[pid][0]
          yvel   = velocity[pid][1]
          zvel   = velocity[pid][2]
          pidf.append(pid)
	  xvelf.append(xvel)
	  yvelf.append(yvel)
	  zvelf.append(zvel)
          #%st = "%d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n"%(pid, xpos, ypos, zpos, xvel, yvel, zvel)
  #typea=(configuration/mpc)%2
  return pidf, xposf, yposf, zposf, xvelf, yvelf, zvelf, box_x, box_y, box_z
  
def runMagic(pidf,  xposf, yposf, zposf, xvelf, yvelf, zvelf, box_x, box_y, box_z,bending,num_chains,monomers_per_chain,nodeGrid,rc,skin,cellGrid,system):
    box=(box_x, box_y, box_z)
    neiListx, neiListy, neiListz = neiListHom(nodeGrid,box,rc,skin)
    
    print 'nei List x',neiListx
    print 'nei List y',neiListy
    print 'nei List z',neiListz
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid,neiListx,neiListy,neiListz)
    mass     = 1.0  
    props    = ['id', 'type', 'mass', 'pos', 'v']
    bondlist = espressopp.FixedPairList(system.storage)
    if bending:
        anglelist = espressopp.FixedTripleList(system.storage)
    for i in range(num_chains):
        	chain = []
        	bonds = []
        	angles = []
        	chtype = 0   #all chains of the same type - thermostat is using id to distinguish hot and cold
        	for k in range(monomers_per_chain):
        		idx  =  i * monomers_per_chain + k
			#print idx
        		part = [ pidf[idx],  chtype, mass, espressopp.Real3D(xposf[idx],yposf[idx],zposf[idx]), espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx])]
        		chain.append(part)
        		if k>0:
        			bonds.append((pidf[idx-1], pidf[idx]))
        		if bending and k>1:
        			angles.append((pidf[idx-2], pidf[idx-1], pidf[idx]))
        	system.storage.addParticles(chain, *props)
        	system.storage.decompose()
        	bondlist.addBonds(bonds)
        	if bending:
        		anglelist.addTriples(angles)
    return system.storage
      
