#!/usr/bin/env python                                                               
# -*- coding: iso-8859-1 -*-                                                        

###########################################################################
#                                                                         #
#  ESPResSo++ Benchmark Python script for a Lennard Jones System          #
#                                                                         #
###########################################################################

import time
import espresso

nsteps      = 10
isteps      = 100
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.005

# set temperature to None for NVE-simulations
temperature = 1.0

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print espresso.Version().info()
print 'Setting up simulation ...'
bonds, angles, x, y, z, Lx, Ly, Lz = espresso.tools.convert.lammps.read('polymer_melt.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = espresso.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=4, ydim=4, zdim=2)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
#system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)
system         = espresso.System()
system.rng     = espresso.esutil.RNG()
system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid=espresso.Int3D(8,5,4)
cellGrid=espresso.Int3D(11,18,11)
neiListx=[0,11,22,33,44,55,66,77,88]
neiListy=[0,18,36,54,72,90] # worked also with [0,23,23] for 2 cores
neiListz=[0,11,22,33,44] # worked also with [0,23,23] for 2 cores
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid, neiListx, neiListy, neiListz)
print "nodeGrid: ",nodeGrid, " cellGrid: ",cellGrid, " NeiList: ",neiListx
integrator     = espresso.integrator.VelocityVerlet(system)  
integrator.dt  = timestep
if (temperature != None):
  thermostat             = espresso.integrator.LangevinThermostat(system)
  thermostat.gamma       = 1.0
  thermostat.temperature = temperature
  integrator.addExtension(thermostat)

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, espresso.Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# Lennard-Jones with Verlet list
vl      = espresso.VerletList(system, cutoff = rc)
potLJ   = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espresso.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espresso.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espresso.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espresso.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)

espresso.tools.analyse.info(system, integrator)
start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  espresso.tools.analyse.info(system, integrator)
end_time = time.clock()
espresso.tools.analyse.info(system, integrator)
espresso.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

