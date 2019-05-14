#  Copyright (C) 2012,2013,2016(H)
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


"""
******************************************
**decomp.py** - Auxiliary python functions
******************************************


*  `nodeGrid(n)`:

    It determines how the processors are distributed and how the cells are arranged.
    `n` - number of processes 

*  `cellGrid(box_size, node_grid, rc, skin)`:

    It returns an appropriate grid of cells.
    
*  `tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.2, precision=0.001)`:

    It tunes the skin size for the current system
    
*  `printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.01)`:
    
    It prints time of running versus skin size in the range [minSkin, maxSkin] with
    the step skinStep
"""


import sys
import espressopp

from espressopp import Int3D
from espressopp.Exceptions import Error

import math
import time
from loadbal import qbicity, changeIndex, halfDecomp, addHsymmetry, adaptNeiList,reDistCellsHom   #,findNodesMS

#def nodeGridAdress(box_size,rc,skin,n,eh_size,ratioMS):
#    totNodesCG,totNodesEH=findNodesMS(node_gridX,totCellsEH,totCellsCG,ratioMS)    # Works only in 1D

# need to factorize n to find optimum nodeGrid
def nodeGrid(box_size,rc,skin,n,eh_size,ratioMS=0):
     # sphereAdr=False is not being implemented
     # adrCenter=[adrCenter[0],adrCenter[1],adrCenter[2]]
     # cg_sizeL=adrCenter[0]-eh_size
     # cg_sizeR=adrCenter[0]+eh_size
     ijkmax = 3*n*n+1
     boxList=[box_size[0],box_size[1],box_size[2]]
     #eh_size=[box_size[0],box_size[1],box_size[2]]   #IDEA or ratioMS[0,1,2]
     if ratioMS>1:
         boxList=[box_size[0]+(2.*eh_size)*ratioMS,box_size[1]*ratioMS,box_size[2]*ratioMS]
     else:
         print "Non AdResS Nodes DD!"
     LoN_Avgmin=sum(boxList)
     ima=boxList.index(max(boxList))
     imi=boxList.index(min(boxList))
     dN = [1,1,1]
     fdN=[0,0,0]
     for i in range(1,n+1):
         for j in range(i,n+1):
             for k in range(j,n+1):
                 if (i*j*k == n) and (i*i + j*j + k*k < ijkmax):
                     dN[0] = k
                     dN[1] = j
                     dN[2] = i
                     ijkmax = i*i + j*j + k*k
                     if qbicity(box_size,rc,skin)==False:
                         ndN=changeIndex(dN,ima,imi)[:]
                         # Adding weighted Averages
                         LoN_norm=[boxList[0]/ndN[0],boxList[1]/ndN[1],boxList[2]/ndN[2]]
                         LoN_Avg=sum(LoN_norm)/3.0
                         if LoN_Avg<=LoN_Avgmin:
                             LoN_Avgmin=LoN_Avg
                             fdN=ndN[:]
                             ijkmax = fdN[0]*fdN[0]+fdN[1]*fdN[1]+fdN[2]*fdN[2]
                             print fdN
                         else:
                             ijkmax = fdN[0]*fdN[0]+fdN[1]*fdN[1]+fdN[2]*fdN[2]
                             print 'No update of dN'
                     else:
                         print 'qbicity check passed DD-UJ'
                         fdN=[k,j,i]
     if abs(box_size[1]-box_size[2])<(2*(rc+skin)):
         if fdN[2]>fdN[1]:
             aux=fdN[2]
             fdN[2]=fdN[1]
             fdN[1]=aux
         else:
             print 'Size Lenghts are eq!'
     else:
         print 'Size Lenghts are different!'
     return Int3D(fdN[0],fdN[1],fdN[2])

def cellGrid(box_size, node_grid, rc,skin):
  rc_skin = rc + skin
  if rc_skin == 0:
    raise Error("interaction range (cutoff + skin) must be larger than 0")
  if (node_grid[0]<=0 or node_grid[1]<=0 or node_grid[2]<=0):
    raise Error("invalid node grid %s" % str(node_grid))
  ix = box_size[0] / (rc_skin * node_grid[0])

  if ix < 1:
    raise Error("local box size in direction 0 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (ix, rc_skin))
  iy = box_size[1] / (rc_skin * node_grid[1])
  if iy < 1:
    raise Error("local box size in direction 1 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (iy, rc_skin))
  iz = box_size[2] / (rc_skin * node_grid[2])
  if iz < 1:
    raise Error("local box size in direction 2 (=%6f) is smaller than interaction range (cutoff + skin = %6f).\n \
                 hint: number of CPUs maybe too high or is prime." % (iz, rc_skin))
  
  return Int3D(ix, iy, iz)

def neiListHom(node_grid,box,rc,skin):
    print "Current NodeGrid",node_grid
    rc_skin=rc+skin
    neiListxin=[]
    neiListyin=[]
    neiListzin=[]
    cursor=[box[0],box[1],box[2]]
    neiListx=[0]*(node_grid[0]+1)
    neiListy=[0]*(node_grid[1]+1)
    neiListz=[0]*(node_grid[2]+1)
    neiListxin=reDistCellsHom(node_grid[0],cursor[0],rc_skin)
    neiListx=adaptNeiList(neiListxin)
    neiListyin=reDistCellsHom(node_grid[1],cursor[1],rc_skin)
    neiListy=adaptNeiList(neiListyin)
    neiListzin=reDistCellsHom(node_grid[2],cursor[2],rc_skin)
    neiListz=adaptNeiList(neiListzin)

    return map(int,neiListx),map(int,neiListy),map(int,neiListz)

def neiListAdress(node_grid, cell_grid,rc,skin,eh_size,adrCenter,ratioMS,idealGasFlag=True,sphereAdr=False):
    # dataStructure Initialization
    print "Current NodeGrid",node_grid
    rc_skin=rc+skin
    #define Neighbor vectors
    neiListx=[0]*(node_grid[0]+1)
    neiListy=[0]*(node_grid[1]+1)
    neiListz=[0]*(node_grid[2]+1)
    # Check adress regions
    adrCenter=[adrCenter[0],adrCenter[1],adrCenter[2]]
    # slab
    cursor=[adrCenter[0]*2,adrCenter[1]*2,adrCenter[2]*2]
    print "cursor",cursor
    cg_sizeR=adrCenter[0]+eh_size
    cg_sizeL=adrCenter[0]-eh_size
    cellsX=round(cg_sizeL/rc_skin-0.5)+round((cursor[0]-cg_sizeR)/ rc_skin-0.5)+round((cg_sizeR-cg_sizeL)/rc_skin-0.5)
    cellsY=round(cursor[1] / rc_skin-0.5)
    cellsZ=round(cursor[2] / rc_skin-0.5)
    print "cells NEWS",cellsX,cellsY,cellsZ
    neiListxin=[]
    neiListyin=[]
    neiListzin=[]
    halfneilListx=[]
    #                               new idea IMPLEMENTATION
    if not sphereAdr:
        #print "All IN:",adrCenter[0],rc_skin,eh_size,int(round(node_grid[0]/2.-0.5))
        halfneilListx=halfDecomp(adrCenter[0],rc_skin,eh_size,int(round(node_grid[0]/2.-0.5)),True)
        print "My halfneilListx...,",halfneilListx
        neiListxin=addHsymmetry(halfneilListx,eh_size,rc_skin,node_grid[0],cellsX,ratioMS,True)  #addHsymmetry(halfNeilListX,eh_size,rc_skin,node_gridX,cellsX,ratioMS)
        print "My neiListxin...",neiListxin         
        neiListx=adaptNeiList(neiListxin)
        ##### OLD but worked
        neiListyin=reDistCellsHom(node_grid[1],cursor[1],rc_skin)
        neiListy=adaptNeiList(neiListyin)
        neiListzin=reDistCellsHom(node_grid[2],cursor[2],rc_skin)
        neiListz=adaptNeiList(neiListzin)
    elif sphereAdr:
        print "WARNING: spherical symmetry is not yet implemented...nevertheless it does not represent a big deal!"
        #quarterneilListx=halfDecomp(adrCenter[0],rc_skin,eh_size,round(node_grid[0]/4-0.5),idealGas)
        #quarterneilListx=halfDecomp(adrCenter[0],rc_skin,eh_size,round(node_grid[0]/4-0.5),idealGas)
    return Int3D(cellsX,cellsY,cellsZ),map(int,neiListx),map(int,neiListy),map(int,neiListz)

######+++++++++++++++++++CONTINUE HERE


'''
    if (idealGasFlag):
        i=1
        f=1
        #exCells=round((cg_sizeR-cg_sizeL)/ rc_skin-0.5)-2
        #jjj=1
        for i in range(1,node_grid[0]+1):
            if (cursor[0] > cg_sizeR):#round(cg_sizeL/ rc_skin-0.5)-cellsX < 0:
                print "ehsize",eh_size
                neiListx[i]=round(cg_sizeL/ rc_skin-0.5)+1
                cursor[0]=cursor[0]-cg_sizeL#-2*rc_skin #test
                neiListx[i]=neiListx[i]+neiListx[i-1]
                print "My i",i
                print"KURSOR ",cursor[0]
                #node_grid[0]=node_grid[0]-1
            elif (cursor[0] >= cg_sizeL) and (cursor[0] <= cg_sizeR) and cursor[0]>0: #(round(adrCenter[0]/ rc_skin-0.5)-cellsX < 0) and (round((cg_sizeR-cg_sizeL)/ rc_skin-0.5)>0):
                #nodeGrid(cg_sizeR-cg_sizeL,node_grid,rc,skin)
                #neRef=node_grid[0]
                deltaNodes=(round((cg_sizeR-cg_sizeL)/rc_skin-0.5)-2)-node_grid[0]-1
                if ((round((cg_sizeR-cg_sizeL)/rc_skin-0.5)-2)%(node_grid[0]-1)!=0):
                    for jjj in range(i,node_grid[0]):
                        #print "fucking i",i
                        neiListx[jjj]=1+neiListx[jjj-1]
                        #cursor[0]=cursor[0]-1*rc_skin
                        if jjj==i:
                            neiListx[jjj]=neiListx[jjj]+1
                            #cursor[0]=cursor[0]-(2)*rc_skin
                        elif jjj==round(node_grid[0]/2):
                            neiListx[jjj]=neiListx[jjj]+1#deltaNodes/deltaNodes
                        elif jjj==round(node_grid[0]/2+1):
                            neiListx[jjj]=neiListx[jjj]+1
                            #cursor[0]=cursor[0]-(2)*rc_skin
                        elif jjj==node_grid[0]-1:
                            neiListx[jjj]=neiListx[jjj]+1#deltaNodes/deltaNodes
                            #cursor[0]=cursor[0]-(2)*rc_skin
                        #else:

                            #neiListx[j]=1+neiListx[j-1]
                        #i=i+1
                        deltaNodes=deltaNodes-1
                    cursor[0]=cursor[0]-(cg_sizeR-cg_sizeL)-rc_skin
                    print "cursor",cursor[0]#=cursor[0]-neiListx[j]*rc_skin
                    print neiListx
                    #i=jjj
                    #cursor[0]=cursor[0]-(1)*rc_skin # extra
                    #print
                    #if deltaNodes<0:
                    #    deltaNodes=0
                    #neiListx[i]=1+deltaNodes+neiListx[i-1]
                    #deltaNodes=deltaNodes-1
                    #cursor[0]=cursor[0]-(1+deltaNodes)*rc_skin #test
                else: #whatever Not implemented
                    print "BADAASSS"
                    #for j in range(node_grid[0]-i,1,-1):
                    #    neiListx[i]=round((cursor[0]-cg_sizeL)/ rc_skin*j-0.5)# change /j  # Homogeneously distributed
                    #    neiListx[i]=neiListx[i]+neiListx[i-1]
                    #    cursor[0]=cursor[0]-(cursor[0]-cg_sizeL)/j
                    #    i=i+1
                #node_grid[0]=node_grid[0]-len(range(node_grid[0],1,-1))

            elif ((cursor[0] < cg_sizeL) and cursor[0]>0):#round(cg_sizeL/ rc_skin-0.5)-cellsX < 0:
                print "ehsize",eh_size
                neiListx[jjj+1]=round(cg_sizeL/ rc_skin-0.5)+1
                cursor[0]=cursor[0]-cg_sizeL#-2*rc_skin #test
                neiListx[jjj+1]=neiListx[jjj+1]+neiListx[jjj]
                print "My i",i
                print"KURSOR ",cursor[0]#i=jjj
    elif f==0:
        #To be implemented soon
        for ik in range(1,node_grid[0]+1):
            if (cursor[0] < cg_sizeL) or (cursor[0] > cg_sizeR):#round(cg_sizeL/ rc_skin-0.5)-cellsX < 0:
                neiListx[ik]=round(cg_sizeL/ rc_skin-0.5)+1
                cursor[0]=cursor[0]-cg_sizeL
                neiListx[ik]=neiListx[ik]+neiListx[i-1]
                node_grid[0]=node_grid[0]-1
            elif (cursor[0] > cg_sizeL) and (cursor[0] < cg_sizeR): #(round(adrCenter[0]/ rc_skin-0.5)-cellsX < 0) and (round((cg_sizeR-cg_sizeL)/ rc_skin-0.5)>0):
                #nodeGrid(cg_sizeR-cg_sizeL,node_grid,rc,skin)
                ka=ik
                for j in range(node_grid[0],1,-1):
                    neiListx[ka]=round((cursor[0]-cg_sizeL)/ rc_skin-0.5)/j  # Homogeneously distributed
                    neiListx[ka]=neiListx[ka]+neiListx[ka-1]-1
                    cursor[0]=cursor[0]-(cursor[0]-cg_sizeL)/j
                    ka=ka+1
                node_grid[0]=node_grid[0]-len(range(node_grid[0],1,-1))
    ia=1
    for ii in range(node_grid[1],0,-1):
        neiListy[ia]=round((cursor[1])/ rc_skin-0.5)/ii  # Homogeneously distributed
        neiListy[ia]=neiListy[ia]+neiListy[ia-1]
        cursor[1]=cursor[1]-(cursor[1])/ii
        ia=ia+1

    ja=1
    for jj in range(node_grid[2],0,-1):
        neiListz[ja]=round((cursor[2])/ rc_skin-0.5)/jj  # Homogeneously distributed
        neiListz[ja]=neiListz[ja]+neiListz[ja-1]
        cursor[2]=cursor[2]-(cursor[2])/jj
        ja=ja+1


    return Int3D(cellsX,cellsY,cellsZ),map(int,neiListx),map(int,neiListy),map(int,neiListz)
'''

def tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.5, precision=0.001, printInfo=True):
  if printInfo:
    print 'The tuning is started. It can take some time depending on your system.'
  
  fi = (1.0+math.sqrt(5.0))/2.0 # golden ratio
  
  npart = espressopp.analysis.NPart(system).compute()
  
  # this is an empirical formula in order to get the appropriate number of steps
  nsteps = int( espressopp.MPI.COMM_WORLD.size * 1000000.0 / float(npart) )
  
  if printInfo:
    print 'CellGrid before tuning: ', system.storage.getCellGrid()
    sys.stdout.write('\nSteps     = %d\n' % nsteps)
    sys.stdout.write('Precision = %g\n' % precision)
    sys.stdout.write('It runs till deltaSkin<precision\n')
  
  if printInfo:
    prnt_format1 = '\n%9s %10s %10s %10s %14s\n'
    sys.stdout.write(prnt_format1 % ('time1: ',' time2: ',' skin1: ', ' skin2: ', ' deltaSkin: '))
  
  while (maxSkin-minSkin>=precision):
    skin1 = maxSkin - (maxSkin-minSkin)/fi
    skin2 = minSkin + (maxSkin-minSkin)/fi

    system.skin = skin1
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time1 = end_time - start_time

    system.skin = skin2
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time2 = end_time - start_time

    if(time1>time2):
      minSkin = skin1
    else:
      maxSkin = skin2
      
    if printInfo:
      prnt_format2 = '%7.3f %10.3f %11.4f %10.4f %12.6f\n'
      sys.stdout.write(prnt_format2 % (time1, time2, minSkin, maxSkin, (maxSkin-minSkin)) )
    
      sys.stdout.write('\nNew skin: %g\n' % system.skin)
      sys.stdout.write('\nNew cell grid: %s\n' % system.storage.getCellGrid())
  
  system.skin = (maxSkin+minSkin)/2.0
  system.storage.cellAdjust()
  
  return (maxSkin+minSkin)/2.0

def printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.005):
  npart = espressopp.analysis.NPart(system).compute()
  # this is an empirical formula in order to get the appropriate number of steps
  nsteps = int( espressopp.MPI.COMM_WORLD.size * 20000000.0 / float(npart) )
  
  print '      Calculations is started. It will print out the dependece of time of \n\
      running of %d steps on the skin size into the file \'timeVSskin.dat\'.\n\
      The range of skin sizes is [%g, %g], skin step is %g. It can take some \n\
      time depending on your system.' % (nsteps, minSkin, maxSkin, skinStep)
  
  curSkin = minSkin
  
  fmt2 = ' %8.4f %8.4f\n' # format for writing the data
  nameFile = 'timeVSskin.dat'
  resFile = open (nameFile, 'w')

  count = 0
  while (curSkin < maxSkin):
    system.skin = curSkin
    system.storage.cellAdjust()
    start_time = time.time()
    integrator.run(nsteps)
    end_time = time.time()
    time1 = end_time - start_time
    
    resFile.write(fmt2 % ( system.skin, time1 ))

    count = count +1
    if (count == 20):
      print 'skin: ', system.skin
      count = 0
    
    curSkin = curSkin + skinStep
  
  resFile.close()
  
  return
  
