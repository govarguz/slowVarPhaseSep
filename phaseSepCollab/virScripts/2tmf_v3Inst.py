#test to load melt data
import espressopp
#import "_espressopp.so"
import mpi4py.MPI as MPI
import sys
import os
import time
from espressopp.tools.decomp import neiListHom
import csv
from magic import runMagic,myGatherConfOTF

num_chains = int(sys.argv[1])
monomers_per_chain = int(sys.argv[2])
temperature1 = 1.0
temperature2 = float(sys.argv[3])
rc=pow(2.0,1.0/6.0)
skin=0.5
dt = float(sys.argv[4])
bending = int(sys.argv[5])
tin = int(sys.argv[6])

if bending:
	kname = 'k1p5'
	ktheta = 1.5# /ptmp/gvargas/work/collabJSmrek/data/
        inconfpath = '/data/isilon/vargas/phaseSep/e++ShWallsInst201606/phaseSepCollab/configurations/N'+str(monomers_per_chain)+'/'
        inconffile = inconfpath+'prc_nc'+str(num_chains)+'N'+str(monomers_per_chain)+kname+'dT'+str(temperature2)+'t'+str(tin)+'.dat'
	outconfpath = '/data/isilon/vargas/phaseSep/e++ShWallsInst201606/phaseSepCollab/configurations/N'+str(monomers_per_chain)+'/out'

'''
else:
	kname = 'kzero'
	ktheta = 0.0
	inconfpath = '/ptmp/jsmrek/Work/phaseg/simulations/data/configurations/N'+str(monomers_per_chain)+'/'
	inconffile = inconfpath+'prc_nc'+str(num_chains)+'N'+str(monomers_per_chain)+kname+'dT'+str(temperature2)+'t'+str(tin)+'.dat'
	outconfpath = '/ptmp/jsmrek/Work/phaseg/simulations/data/configurations/N'+str(monomers_per_chain)+'/'#'/ptmp/jsmrek/Work/phaseg/simulations/runscripts/local_test/'

'''

os.system('mkdir '+outconfpath)
#### get here!

pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espressopp.tools.readxyz(inconffile)
box = (Lxf, Lyf, Lzf)
print "Box Size is:",box
system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(box,rc,skin,MPI.COMM_WORLD.size,0)
print "your NodeGrid is",nodeGrid
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)

neiListx, neiListy, neiListz = neiListHom(nodeGrid,box,rc,skin)

print 'nei List x',neiListx
print 'nei List y',neiListy
print 'nei List z',neiListz

system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid,neiListx,neiListy,neiListz)
  
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


fid = 'nc'+str(num_chains)+'N'+str(monomers_per_chain)+kname+'dT'+str(temperature2)

#PRINTING PARAMETERS TO CHECK
print 'PARAMETERS OF THE RUN'
print '---------------------'
os.system('hostname')
print '{:22}{:25}'.format('infile:', inconffile)
print '{:22}{:25}'.format('num_chains:', num_chains)
print '{:22}{:25}'.format('monomers per chain:', monomers_per_chain)
print '{:22}{:25}'.format('temperature 2:', temperature2)
print '{:22}{:25}'.format('k_theta:', ktheta)
print '{:22}{:25}{:30}'.format('thermostat gammas:', thermostat.gamma1, thermostat.gamma2)
print '{:22}{:25}'.format('time step:', dt)
print '{:22}{:25}'.format('fid:', fid)
print '{:22}{:25}'.format('outpath:', outconfpath)

#SETTING UP CONFORMATION
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

#INTERACTIONS SETUP
#same LJ interaction between all particles
#LJ nonbonded
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


intime = time.time()
tsamp = 10000 #tau
t=tin

# H's hack
f1=open("NIntPerRank.dlb","wt")
wr1=csv.writer(f1,delimiter=" ")
pt1=[]
f2=open("NPartsPerRank.dlb","wt")
wr2=csv.writer(f2,delimiter=" ")
pt2=[]
f3=open("PartsHangoverPerRankPerDir.dlb","wt")
wr3=csv.writer(f3,delimiter=" ")
pt3=[]
f4=open("PartsCommCellsPerRankPerDir.dlb","wt")
wr4=csv.writer(f4,delimiter=" ")
pt4=[]
'''
Finished creating load balancing Ranks
'''
dd=True
espressopp.tools.myfastwritexyz(outconfpath+'prc_'+fid+'t'+str(int(t))+'.dat',system, velocities=True, unfolded = True,mpc = monomers_per_chain)# mywrite writes the type as particleId/monomers_per_chain %2
nsteps = 3#200*int(thermostat.gamma1) #just a large number
tausteps = int(1.0/dt)
mdSteps=0

for j in range(nsteps): 		#run nsteps x samp time
    espressopp.tools.analyse.info(system,integrator)
    integrator.run(10)#tsamp*tausteps) 		#run for samp time * \tau
    if mdSteps>9:#3.9999e6: #Tunned Param for J Smrek
        l=len(system.getInteraction(0).getVerletList().getAllPairs()) 
        p=[]
        for k2 in range(l):
            p.append(len(system.getInteraction(0).getVerletList().getAllPairs()[k2]))
        pt1.append(p)
        pt2.append(espressopp.analysis.NRealPart(system).compute())
        pt3.append(system.storage.getCounters())
        pt4.append(system.storage.getCommCellsTotalParticles())
        system.storage.resetParticlesHangoverCounters()
        system.storage.resetParticlesCommCellsCounters()
        t = integrator.step*integrator.dt + tin
        mdSteps=mdSteps+tsamp*tausteps
        #espressopp.tools.myfastwritexyz(outconfpath+'prc_'+fid+'t'+str(int(t))+'.dat',system, velocities=True, unfolded = True, mpc = monomers_per_chain)
        if dd:
	    xposf2=[]	
            pidf2, xposf2, yposf2, zposf2, xvelf2, yvelf2, zvelf2, box_x, box_y, box_z=myGatherConfOTF(system,velocities=True,unfolded=True,append=False,mpc=monomers_per_chain)
	    #print 'print Xposf',xposf2[2]
            runMagic(pidf2, xposf2, yposf2, zposf2, xvelf2, yvelf2, zvelf2, box_x, box_y, box_z,bending,num_chains,monomers_per_chain,nodeGrid,rc,skin,cellGrid,system)
            print "There you go, Magic has started"
        else:
            print "4M but no additional DD is required!"
    else:
        print "No DD update"
    t = integrator.step*integrator.dt + tin
    mdSteps=mdSteps+tsamp*tausteps
    espressopp.tools.myfastwritexyz(outconfpath+'prc_'+fid+'t'+str(int(t))+'.dat',system, velocities=True, unfolded = True, mpc = monomers_per_chain)

wr1.writerows(pt1)
wr2.writerows(pt2)
wr3.writerows(pt3)
wr4.writerows(pt4)
f1.close()
f2.close()
f3.close()
f4.close()

fintime = time.time()
espressopp.tools.analyse.final_info(system, integrator, vl, intime, fintime)

print "Done HAVG."
print "Running time: ", fintime-intime











