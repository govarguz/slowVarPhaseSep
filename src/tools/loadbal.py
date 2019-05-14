# -*- coding: utf-8 -*-
"""
******************************************
**loadbal.py** - Auxiliary python functions
******************************************


*  `qbicity(box_size,rc,skin)`:
    
    Its a function to check the system size cubicity, with a tolerance given by the rc+skin
    `dLnorm` - (Lx,Ly,Lz)/Lmax 

*  `cellGrid(box_size, node_grid, rc, skin)`:

    It returns an appropriate grid of cells.
    
*  `tuneSkin(system, integrator, minSkin=0.01, maxSkin=1.2, precision=0.001)`:

    It tunes the skin size for the current system
    
*  `printTimeVsSkin(system, integrator, minSkin=0.01, maxSkin=1.5, skinStep = 0.01)`:
    
    It prints time of running versus skin size in the range [minSkin, maxSkin] with
    the step skinStep

@author: vargas
"""
import sys
import espressopp

from espressopp import Int3D, Real3D
from espressopp.Exceptions import Error
from random import randint

def qbicity(box_size,rc,skin,cellDomFactor=2.0):
    rc_skin=rc+skin
    au=[box_size[0],box_size[1],box_size[2]]
    ima=au.index(max(au))
    imi=au.index(min(au))
    if (box_size[ima]-box_size[imi])<(cellDomFactor*rc_skin):
        print 'qbicity check passed'
        qFlag=bool(1)
    else:
        qFlag=bool(0)     
    return qFlag

def changeIndex(dN,ima,imi):
    '''
    ima correspond to the index of the maximum size of the box and imi to the minimum
    one use push
    '''
    print dN
    aux=[0,0,0]
    ndN=[0,0,0]
    dima=dN.index(max(dN))
    dimi=dN.index(min(dN))
    ndN[ima]=dN[dima]
    ndN[imi]=dN[dimi]
    listInd=range(3)
    if ima>imi:
        listInd.pop(ima)
        listInd.pop(imi)
    else:
        listInd.pop(imi)
        listInd.pop(ima)
    aux=dN[:]
    if dima>dimi:  
        aux.pop(dima)
        aux.pop(dimi)
    else:
        aux.pop(dimi)
        aux.pop(dima)
    ndN[listInd[0]]=aux[0]
    return ndN

def halfDecomp(adrCenter1D,rc_skin,eh_size,halfCores1D,idealGas):
    '''
    this function decomposes half of the box (ojo: in 1 direction ONLY)
    '''
    pLoadIG=1
    cellSizes=[]
    usedCores=0
    #flagHD=0
    for i in range(halfCores1D):
        if idealGas:
            if i==0: #and ((halfCores1D-usedCores)<round((eh_size)/rc_skin-0.5)): # Added 2nd condition
                cellSizes.append(round((adrCenter1D-eh_size)/rc_skin-0.5)+pLoadIG)
                usedCores=1
                #flagHD=1
            else:
                print "works"
                #usedCoresFix=usedCores
                #while (halfCores1D-usedCores)>round((eh_size)/rc_skin-0.5):
                #    cellSizes.append(round((adrCenter1D-eh_size)/rc_skin/(halfCores1D-usedCoresFix-round((eh_size)/rc_skin-0.5))-0.5))
                #    usedCores=+1
                [cellSizes.append(round((eh_size)/rc_skin/(halfCores1D-usedCores)-0.5)) for i in range(usedCores,halfCores1D)]
                deltaCells=round((eh_size)/rc_skin-0.5)-round((eh_size)/rc_skin/(halfCores1D-usedCores)-0.5)*(halfCores1D-usedCores)
                #if flagHD:
                #    cellSizes[usedCores]=cellSizes[usedCores]+deltaCells-pLoadIG
                #else:
                #    cellSizes[usedCores]=cellSizes[usedCores]+deltaCells
                #halfCores1D=i
                cellSizes[usedCores]=cellSizes[usedCores]+deltaCells-pLoadIG
                return cellSizes
                #cellSizes.append(round((adrCenter1D-eh_size)/rc_skin-0.5))
        else: # TO BE MOD H's
            if i==0:
                cellSizes.append(round((adrCenter1D-eh_size)/rc_skin-0.5))
                usedCores=1
            else:
                print "NOT YET working"
                [cellSizes.append(round((eh_size)/rc_skin-0.5)) for i in range(usedCores,halfCores1D)]
                halfCores1D=i
    #cellSizes.append(round((adrCenter1D-eh_size)/rc_skin-0.5)+pLoadIG)
    return cellSizes

def addHsymmetry(halfNeilListX,eh_size,rc_skin,node_gridX,cellsX,ratioMS,idealGas):
    '''
    this function decomposes half of the box (ojo: in 1 direction ONLY)
    '''
    #neiListX=[0]*(node_gridX)
    wholeNeilListX=[]
    aux=halfNeilListX[:]
    aux.reverse()
    halfNeilListX.extend(aux)
    aux2=0
    if len(halfNeilListX)<node_gridX:
        if (node_gridX-len(halfNeilListX))==1:
            print "Number of cores is not even!"
            aux2=halfNeilListX[len(halfNeilListX)-1]
            halfNeilListX.append(round(halfNeilListX[len(halfNeilListX)-1]/2.-0.5))
            halfNeilListX[len(halfNeilListX)-2]=aux2-halfNeilListX[len(halfNeilListX)-1]
            if sum(halfNeilListX)!=cellsX:
                print "Error in H's initial DD  HACK"
                #diffCell=cellsX-sum(halfNeilListX)
                wholeNeilListX=reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas) #halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS
            else:
                if any([ v == 0 for v in halfNeilListX]):                # Recently added 138, 139 and 140    
                    wholeNeilListX=reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas)
                else:    
                    print "addHsymmetry all tests passed"
                    wholeNeilListX=halfNeilListX[:]
            #wholeNeilListX=halfNeilListX[:]
        else:
            print "Error in H's initial DD  HACK"
    elif len(halfNeilListX)==node_gridX and aux2==0:
        if sum(halfNeilListX)!=cellsX:
            print "Error in H's initial DD  HACK"
            #diffCell=cellsX-sum(halfNeilListX)
            wholeNeilListX=reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas)
        else:
            if any([ v == 0 for v in halfNeilListX]):                    # Recently added 152, 153 and 154
                wholeNeilListX=reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas)
            else:    
                print "addHsymmetry all tests passed"
                wholeNeilListX=halfNeilListX[:]
    return wholeNeilListX
    
def adaptNeiList(neiListxin):
    neiListx=[]
    neiListx.append(0)
    [neiListx.append(neiListxin[i]+neiListx[i]) for i in range(len(neiListxin)-1)]
    neiListx.append(neiListxin[len(neiListxin)-1]+neiListx[len(neiListx)-1])
    print "Your subscribed Neigh List IS:",neiListx
    return neiListx          

def reDistCellsHom(node_gridX,sizeX,rc_skin):
    wholeNeiListX=[]
    cellsX=round(sizeX/rc_skin-0.5)
    if node_gridX%2==0 and cellsX%2==0:
        [wholeNeiListX.append(cellsX/node_gridX) for i in range(node_gridX)]
    elif node_gridX%2!=0 and cellsX%2!=0:
        [wholeNeiListX.append(round((cellsX)/node_gridX-0.5)) for i in range(node_gridX)]
        if int(cellsX-sum(wholeNeiListX))!=0:
            wholeNeiListX=redistDeltaRandomly(wholeNeiListX,cellsX-sum(wholeNeiListX),0) # passing Delta
        else:
            print "You are Well.. DONE!"            
    else:
        if node_gridX%2==0:
            [wholeNeiListX.append((cellsX-1)/node_gridX) for i in range(node_gridX)]
            wholeNeiListX[node_gridX-1]=wholeNeiListX[node_gridX-1]+1   # Punishing the last one
        elif cellsX%2==0:
            [wholeNeiListX.append(round((cellsX)/node_gridX-0.5)) for i in range(node_gridX)]
            wholeNeiListX=redistDeltaRandomly(wholeNeiListX,cellsX-sum(wholeNeiListX),0) # passing Delta    
    return wholeNeiListX


def reDistCells(halfNeilListX,cellsX,eh_size,rc_skin,node_gridX,ratioMS,idealGas):
    '''This function is matching proportion of cells to the
    cores on a region basis '''
    print "WARNING: This only works whenever the Cells1D are at least twice as big as Nodes1D!!!"
    wholeNeiListX=[]
    totCellsEH=round(2*eh_size/rc_skin-0.5)
    totCellsCG=cellsX-totCellsEH
    #ratioCells=1.0*totCellsEH/(cellsX)
    if ratioMS==1:
        if node_gridX%2==0 and cellsX%2==0:
            [wholeNeiListX.append(cellsX/node_gridX) for i in range(node_gridX)]
        else:
            if node_gridX%2==0:
                [wholeNeiListX.append((cellsX-1)/node_gridX) for i in range(node_gridX)]
                wholeNeiListX[node_gridX-1]=+1
            if cellsX%2==0:
                [wholeNeiListX.append(round((cellsX)/node_gridX-0.5)) for i in range(node_gridX)]
                wholeNeiListX=redistDeltaRandomly(wholeNeiListX,cellsX-sum(wholeNeiListX),0) # passing Delta
    else:
        totNodesCG,totNodesEH=findNodesMS(node_gridX,totCellsEH,totCellsCG,ratioMS,idealGas)  #added all int() functions
        print "NOdes CG and Nodes EH:",totNodesCG,totNodesEH
        if node_gridX%2==0 and cellsX%2==0:
            wholeNeiListX_EH=[0]*(node_gridX)
            wholeNeiListX_CG=[0]*(node_gridX)
            wholeNeiListX=[0]*(node_gridX)
            if totNodesCG%2==0:
                indCG1=int(totNodesCG)/2
                ##### CONTINUE here!!!
                indEH1=indCG1+int(totNodesEH)
                au1=range(int(indCG1))
                au1.extend(range(indEH1,indEH1+indCG1))
                for i in au1:
                    wholeNeiListX_CG[i]=round(1.0*totCellsCG/totNodesCG-0.5)
                print "wholeNeiListX_CG is",wholeNeiListX_CG
                if totNodesEH%2==0:
                    for i in range(indCG1,indEH1):
                        wholeNeiListX_EH[i]=round(1.0*totCellsEH/totNodesEH-0.5)
                else:
                    for i in range(indCG1,indEH1):
                        wholeNeiListX_EH[i]=round(1.0*(totCellsEH-1)/totNodesEH-0.5)
                    wholeNeiListX[indEH1-1]=+1  # Punishing last Node
                for k in range(node_gridX):    
                    wholeNeiListX[k]=wholeNeiListX_EH[k]+wholeNeiListX_CG[k]    
            else:  #TO BE IMPROVED not fullfilling all nodes!!!
                indCG1=int(totNodesCG/2.0)
                indEH1=indCG1+int(totNodesEH)
                au2=range(int(indCG1))
                au2.extend(range(indEH1,indEH1+indCG1))
                for i in au2:
                    wholeNeiListX_CG[i]=round(1.0*(totCellsCG-1)/totNodesCG-0.5)
                wholeNeiListX_CG[indEH1+indCG1-1]=+1  # Punishing last Node
                if totNodesEH%2==0:
                    for i in range(indCG1,indEH1):
                        wholeNeiListX_EH[i]=round(1.0*totCellsEH/totNodesEH-0.5)
                else:
                    for i in range(indCG1,indEH1):
                        wholeNeiListX_EH[i]=round(1.0*(totCellsEH-1)/totNodesEH-0.5)
                    wholeNeiListX_EH[indEH1-1]=+1  # Punishing last Node
                for k in range(node_gridX):    
                    wholeNeiListX[k]=wholeNeiListX_EH[k]+wholeNeiListX_CG[k]    
        else: # TO BE IMPROVED
            if node_gridX%2==0:
                [wholeNeiListX.append(round((cellsX-1)/node_gridX-0.5)) for i in range(node_gridX)]   # Homogeneously distributed
                wholeNeiListX[node_gridX-1]=+1
            if cellsX%2==0:
                [wholeNeiListX.append(round((cellsX)/node_gridX-0.5)) for i in range(node_gridX)]
                wholeNeiListX=redistDeltaRandomly(wholeNeiListX,cellsX-sum(wholeNeiListX),totNodesEH)
    print "My Redist WholeNeiList is:",wholeNeiListX            
    return wholeNeiListX
    
def redistDeltaRandomly(wholeNeiListX,deltaCells,totNodesEH=0):
    wholeNeiListXcopy=wholeNeiListX[:]
    index=len(wholeNeiListX)-1
    indexOut=[0]*int(deltaCells)
    if totNodesEH==0:
        for p in range(0,int(deltaCells)):
            aux2=randint(0,index)
            while randint(0,index)==aux2:
                aux2=randint(0,index)
            indexOut[p]=aux2    
        for i in indexOut:
            wholeNeiListXcopy[i]=wholeNeiListX[i]+1
    else:
        for p in range(0,int(deltaCells)):
            index=len(wholeNeiListX)-1-totNodesEH
            aux2=randint(0,index)
            while randint(0,index)==aux2:
                aux2=randint(0,index)
            indexOut[p]=aux2    
        for i in indexOut:
            wholeNeiListXcopy[i]=wholeNeiListX[i]+1                
    return wholeNeiListXcopy

def findNodesMS(node_gridX,totCellsEH,totCellsCG,ratioMS,idealGas): #Added all any statements and its corresponding IF's    
    fRatioCG=ratioMS*(1.0*totCellsEH/(totCellsCG))
    fRatioEH=(1./ratioMS)*(1.0*totCellsCG/(totCellsCG))
    if idealGas:
        if (node_gridX-2.)<=totCellsEH and node_gridX>totCellsEH:
            totNodesEH=round(totCellsEH/2.-0.5)# Must be tuned for every case!
            totNodesCG=node_gridX-totNodesEH
        else:
            totNodesEH=totCellsEH# Must be tuned for every case!
            totNodesCG=node_gridX-totNodesEH
            print "Not enough nodes..WTE"
            
    else:
        if node_gridX<=(totCellsEH+totCellsCG):
            totNodesEH=round(fRatioEH/(fRatioEH+fRatioCG)*node_gridX)
            totNodesCG=round(fRatioCG/(fRatioEH+fRatioCG)*node_gridX)
            if (totNodesEH+totNodesCG)!=node_gridX:
                if totNodesEH>(totCellsEH):
                    print "NO GO"
                    diffNodesCells=totNodesEH-totCellsEH
                    if diffNodesCells<=(totCellsCG-totNodesCG):
                        print "updating Nodes"
                        totNodesCG=totNodesCG+diffNodesCells                    
                        totNodesEH=totNodesEH-diffNodesCells
                    else:
                        print "Erroneous DD values Nodes more than Cells EH!"
                elif totNodesCG>(totCellsCG):
                    print "NO GO"
                    diffNodesCells=totNodesCG-totCellsCG
                    if diffNodesCells<=(totCellsEH-totNodesEH):
                        print "updating Nodes"
                        totNodesCG=totNodesCG-diffNodesCells
                        totNodesEH=totNodesEH+diffNodesCells
                    else:
                        print "Erroneous DD values Nodes more than Cells CG!"
       
            else:
                print "Everything is fine!"
        else:
            print "Erroneous DD values... Nodes more than Cells!"       
    return totNodesCG,totNodesEH


#
#def chNodeGridGeo(dsum,box_size,ima,midL,cellAvg,qFlag,facCell):
#    facCell=1.1    
#    if qFlag=False and (dsum)>3:
#        box=[box_size[0],box_size[1],box_size[2]]
#        #boxMax=box[box.index(max(box))]
#        #boxNorm=box/(rc+skin)
#    if (box[midL]+1.1*cellAvg)>box[ima]:
#        # iterate again
#    else:
#        # goto RETURN
#        
#     cellLinAvg=int((box_size[0]+box_size[1]+box_size[2])/n/(rc+skin))
