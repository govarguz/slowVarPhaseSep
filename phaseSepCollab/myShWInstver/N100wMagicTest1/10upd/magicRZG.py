# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:46:19 2016

@author: vargas
"""

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
import numpy as np

# Magic is still wired!

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

  #numParticles  = int(espressopp.analysis.NPart(system).compute())
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
  
def runMagic(pidf,  xposf, yposf, zposf, xvelf, yvelf, zvelf, box_x, box_y, box_z,bending,num_chains,monomers_per_chain,nodeGrid,rc,skin,cellGrid,system): #Nikulin
    box=(box_x, box_y, box_z)
    neiListx, neiListy, neiListz = neiListHom(nodeGrid,box,rc,skin)
     
    print 'nei List x',neiListx
    print 'nei List y',neiListy
    print 'nei List z',neiListz

    #espressopp.FixedPairList.remove()
    #espressopp.FixedTripleList.remove()
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
		#integrator     = espressopp.integrator.VelocityVerlet(system)  
		#integrator.dt  = dt
#    return system.storage

def runMagic2(pidf,  xposf, yposf, zposf, xvelf, yvelf, zvelf, box_x, box_y, box_z,bending,num_chains,monomers_per_chain,nodeGrid,rc,skin,cellGrid,ktheta,dt,temperature1,temperature2,pt1,pt2,pt3,pt4,neiListxIn=[], neiListyIn=[], neiListzIn=[]): #
    box=(box_x, box_y, box_z)
    #ktheta=1.5	
    system         = espressopp.System()
    system.rng     = espressopp.esutil.RNG()
    system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin    = skin
    nodeGrid       = espressopp.tools.decomp.nodeGrid(box,rc,skin,MPI.COMM_WORLD.size,0)
    cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)

    if len(neiListxIn)>1:
        print "Cells Neighbors Input List has been given!"
        coRperPlane=1.9 #2.5  # Optimization to be implemented soon
        drive=0 # This Flag makes the KPIs checks only for Part+Ints or if drive=1 all KPIs
        neiListx, neiListy, neiListz =criteriaInst(pt1,pt2,pt3,pt4,neiListxIn,neiListyIn,neiListzIn,coRperPlane,drive,None,threshNParts=26.,threshNInts=270.,threshNHangs=68.,threshNComms=10)
        #criteriaInst(pt1,pt2,pt3,pt4,numNodes,neiListxIn,neiListyIn,neiListzIn,drive=0,filename=None,threshNParts=0,threshNInts=0,threshNHangs=0,threshNComms=0):
    else:    
        print "Cells Neighbors Input List has NOT been given!"
        neiListx, neiListy, neiListz = neiListHom(nodeGrid,box,rc,skin)
            
    print 'nei List x',neiListx
    print 'nei List y',neiListy
    print 'nei List z',neiListz
  
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid,neiListx,neiListy,neiListz)
    #system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)  
    integrator     = espressopp.integrator.VelocityVerlet(system)  
    integrator.dt  = dt
    #LANGEVIN DOUBLE THERMOSTAT
    thermostat             	= espressopp.integrator.LangevinThermostat2TId(system)
    thermostat.gamma1       = 10.0
    thermostat.gamma2       = 10.0
    thermostat.temperature1 = temperature1
    thermostat.temperature2 = temperature2
    thermostat.type1				= 0
    thermostat.type2				= 1
    thermostat.mpc = monomers_per_chain
    integrator.addExtension(thermostat)  
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
		#integrator     = espressopp.integrator.VelocityVerlet(system)  
		#integrator.dt  = dt
    vl 						 = espressopp.VerletList(system, cutoff=rc)
    interLJ00		     = espressopp.interaction.VerletListLennardJones(vl)
    potLJ00			 		 = espressopp.interaction.LennardJones(epsilon = 1.0, sigma = 1.0, cutoff = rc,shift = 'auto')
    interLJ00.setPotential(type1=0, type2=0, potential = potLJ00)
    system.addInteraction(interLJ00)

    # FENE bonds
    potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
    interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
    system.addInteraction(interFENE)

    #bending interation - only if bending switched on
    # Cosine with FixedTriple list
    if bending:
	potCosine = espressopp.interaction.Cosine(K=ktheta, theta0=3.1415926)
	interCosine = espressopp.interaction.FixedTripleListCosine(system, anglelist, potCosine)
	system.addInteraction(interCosine)

    return system,integrator
    
def criteriaInst(pt1,pt2,pt3,pt4,neiListxIn,neiListyIn,neiListzIn,coRperPlane=1.9,drive=0,filename=None,threshNParts=0,threshNInts=0,threshNHangs=0,threshNComms=0):
    #coRperPlane=2.00
    # periodicityVerletL,periodicityObDifuMDsteps,  To be added soon as complementary parameters
    #nodesMagicNum=numNodes*nNodesCrit   # To be further tuned
    if filename==None:
        pt1C=sum(pt1[i]>threshNInts for i in range(len(pt1)))
        pt2C=sum(pt2[i]>threshNParts for i in range(len(pt2)))
        #pt1C=sum(pt1[0][i]>threshNInts for i in range(len(pt1[0]))) #Thisworks if pt1Delta is not the case
        #pt2C=sum(pt2[0][i]>threshNParts for i in range(len(pt2[0])))
        #pt1C=sum(pt1[i]>threshNInts for i in range(len(pt1))) wrong due to [[]]
        #pt2C=sum(pt2[i]>threshNParts for i in range(len(pt2)))
        if drive==1:
            for j in [0,1,2,3,4,5]:# xDir[0,3] ; yDir[1,4] ; zDir[2,5]
                if j==0:
                    pt3CxL=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))  #check if adding [0] is needed
                    pt4CxL=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))
                elif j==1:
                    pt3CyL=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))            
                    pt4CyL=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))                            
                elif j==2:
                    pt3CzL=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))                        
                    pt4CzL=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))                                        
                elif j==3:
                    pt3CxR=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))                                   
                    pt4CxR=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))                                                   
                elif j==4:
                    pt3CyR=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))                                                
                    pt4CyR=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))                                                                
                else:
                    pt3CzR=sum(pt3[i][j]>threshNHangs for i in range(len(pt3)))                                                        
                    pt4CzR=sum(pt4[i][j]>threshNHangs for i in range(len(pt4)))                                                                                        
        else:            
            print "Only Parts AND/OR Ints driven... Check Nodes loading"
            pt3CxL=0
            pt3CxR=0
            pt3CyL=0
            pt3CyR=0
            pt3CzL=0
            pt3CzR=0
            pt4CxL=0
            pt4CxR=0
            pt4CyL=0
            pt4CyR=0
            pt4CzL=0
            pt4CzR=0
        #drive=0            
        if pt1C>coRperPlane*3.:
            ddFlag=1
            print "WARNING: not implemented yet...(add ifs and control infra...) threshNInts"        
        elif pt2C>coRperPlane*3.: #CHANGE
            ddFlag=1
            pt23D=buildMatrix(len(neiListxIn),len(neiListyIn),len(neiListzIn),pt2)
            #pt23D=buildMatrix(len(neiListxIn),len(neiListyIn),len(neiListzIn),pt2[0])
            neiListx,neiListy,neiListz=callRedDistMatrixNParts(pt23D,neiListxIn,neiListyIn,neiListzIn,coRperPlane,threshNParts)
            #else:
        elif any([pt3CxL,pt3CxR,pt3CyL,pt3CyR,pt3CzL,pt3CzR]):
            print "WARNING: not implemented yet...(add ifs and control infra...) threshNHangs"
            ddFlag=1
            #   else:
        elif any([pt4CxL,pt4CxR,pt4CyL,pt4CyR,pt4CzL,pt4CzR]):
            print "WARNING: not implemented yet...(add ifs and control infra...) threshNCommss"
            ddFlag=1
        else:
            print "No DD needed  ....TBD"
    else:
        print "WARNING: not implemented yet...(add ifs and control infra...) READ Files, etc, etc"
        #threshNParts,threshNInts,threshNHangs,threshNComms=readCriteria(filename)
        if pt1>threshNInts:
            ddFlag=1
        else:
            if pt2>threshNParts:
                ddFlag=1
            else:
                if pt3>threshNHangs:
                    ddFlag=1
                else:
                    if pt4>threshNComms:
                        ddFlag=1
                    else:
                        print "No DD needed, your falg is: ",ddFlag
    return neiListx,neiListy,neiListz

def checkDDdyn(pt1,pt2,threshNInts,threshNParts,coRperPlane):
    # periodicityVerletL,periodicityObDifuMDsteps,  To be added soon as complementary parameters
    #nodesMagicNum=numNodes*nNodesCrit   # To be further tune
    ddFlag=0
    print "Tonces q es pt1..",pt1
    #pt1C=sum(pt1[0][i]>threshNInts for i in range(len(pt1[0])))
    pt1C=sum(pt1[i]>threshNInts for i in range(len(pt1)))
    print "pt1C counter is",pt1C
    #pt2C=sum(pt2[0][i]>threshNParts for i in range(len(pt2[0])))
    pt2C=sum(pt2[i]>threshNParts for i in range(len(pt2)))
    print "pt2C counter is",pt2C
    if pt2C>coRperPlane*3.:
        ddFlag=1
        #pt23D=buildMatrix(len(neiListxIn),len(neiListyIn),len(neiListzIn),pt2)
        #neiListx,neiListy,neiListz=callRedDistMatrixNParts(pt23D,neiListxIn,neiListyIn,neiListzIn,coRperPlane)
    elif pt1C>coRperPlane*3.:
        ddFlag=1
    print "from checkDDdyn DD Flag is:",ddFlag    
    return ddFlag
                        
def callRedDistMatrixNParts(pt23D,neiListxIn,neiListyIn,neiListzIn,coRperPlane,threshNParts):
    neiListx=[0]*(len(neiListxIn))
    neiListy=[0]*(len(neiListyIn))
    neiListz=[0]*(len(neiListzIn))    
    #coR=2.5
    iniNLists=1 # This could be further Tuned or changed TBI
    neiListzC=countPlanes(pt23D,neiListxIn,neiListyIn,neiListzIn,threshNParts) # pt23D plane XY -> Defines neiListZ
    print "neiListzC is:",neiListzC
    neiListyC=countPlanes(pt23D,neiListxIn,neiListzIn,neiListyIn,threshNParts) # pt23D plane XZ -> Defines neiListY
    print "neiListyC is:",neiListyC    
    #neiListxC=countPlanesNonQ(pt23D,neiListyIn,neiListzIn,neiListxIn,threshNParts) # pt23D plane XZ -> Defines neiListX
    neiListxC=countPlanes(pt23D,neiListyIn,neiListzIn,neiListxIn,threshNParts) # pt23D plane XZ -> Defines neiListX
    print "neiListxC is:",neiListxC    
    for i in range(iniNLists,len(neiListxIn)-1): # CHANGE added -1
        if neiListxC[i-1]>((len(neiListxIn)-1)/coRperPlane): #added -1
            neiListx[i]=neiListxIn[i]+1
            if i+1==(len(neiListxIn)-1):
                print "ya Fue"
                neiListx[i+1]=neiListxIn[i+1]
            else:    
                neiListxIn[i+1]=neiListx[i+1]-1
        else:
            neiListx[i]=neiListxIn[i]     # TO be Implemented Fro this case!!!       
            if i+1==(len(neiListxIn)-1):            
                neiListx[i+1]=neiListxIn[i+1]   # Actual Testing
            else:
                print "QTDPC"                                
    for j in range(iniNLists,len(neiListyIn)-1):
        if neiListyC[j-1]>((len(neiListyIn)-1)/coRperPlane):
            neiListy[j]=neiListyIn[j]+1
            if j+1==(len(neiListyIn)-1):
                print "ya Fue"
                neiListy[j+1]=neiListyIn[j+1]                
            else:
                neiListyIn[j+1]=neiListy[j+1]-1
        else:
            neiListx[j]=neiListxIn[j]   
            if j+1==(len(neiListyIn)-1):                    
                neiListy[j+1]=neiListyIn[j+1]                
            else:
                print "QTDPC"                
    for k in range(iniNLists,len(neiListzIn)-1):
        if neiListzC[k-1]>((len(neiListzIn)-1)/coRperPlane):
            neiListz[k]=neiListzIn[k]+1
            if k+1==(len(neiListzIn)-1):
                print "ya Fue"
                neiListz[k+1]=neiListzIn[k+1]                
            else:
                neiListzIn[k+1]=neiListz[k+1]-1
        else:
            neiListx[k]=neiListxIn[k]
            if k+1==(len(neiListzIn)-1):                        
                neiListz[k+1]=neiListzIn[k+1]                
            else:
                print "QTDPC"
    return neiListx,neiListy,neiListz   

def countPlanes(pt23D,neiListxIn,neiListyIn,neiListzIn,threshNParts):
    neiListFlex=[]
    for k in range(len(neiListzIn)-1):
        cter=0
        for j in range(len(neiListyIn)-1):
            for i in range(len(neiListxIn)-1):
                if pt23D[i][j][k]>threshNParts:
                    cter=cter+1
        neiListFlex.append(cter)
    return neiListFlex

def countPlanesNonQ(pt23D,neiListxIn,neiListyIn,neiListzIn,threshNParts):
    neiListFlex=[]
    for k in range(len(neiListzIn)-1):
        cter=0
        for j in range(len(neiListyIn)-1):
            for i in range(len(neiListxIn)-1):
                if pt23D[i][j][k]>threshNParts:
                    cter=cter+1
        neiListFlex.append(cter)
    return neiListFlex
                
def mapIndexToPosition(nodeGridx,nodeGridy,nodeGridz,indexRank):    
    px = indexRank % nodeGridx
    indexRank = indexRank/nodeGridx
    py = indexRank % nodeGridy
    indexRank = indexRank/nodeGridy
    pz= indexRank
    return px,py,pz

def buildMatrix(neiListxSize,neiListySize,neiListzSize,matProp): # matProp Is the list of arrays | indexRank  
    indexRank=len(matProp)
    nodeGridx=neiListxSize-1
    nodeGridy=neiListySize-1
    nodeGridz=neiListzSize-1
    ranksMat=np.zeros((nodeGridx,nodeGridy,nodeGridz))
    for i in range(indexRank):
        px,py,pz=mapIndexToPosition(nodeGridx,nodeGridy,nodeGridz,i)
        ranksMat[px,py,pz]=matProp[i]
    return ranksMat