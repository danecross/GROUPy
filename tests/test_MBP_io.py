
from Harvest import *
from Halo import *
from Particle import *
from Energies import *

out_path = "test_data/"
rs = [out_path+"out_%i.list"%i for i in range(10)]
particles = [out_path+"halos_snapshot_%i.hdf5.0.particles"%i for i in range(10)]

#################################
## populate particle directory ##
#################################

pdir = "test_output/particles/"
if not os.path.exists(pdir): os.mkdir(pdir)

harvest_particles(particles[4], 4, 10, pdir)


##############################
## test get particle masses ##
##############################

particle_mass = get_particle_mass(rs[0])
assert(particle_mass==4.61528e+06)

#########################
## test get_halo_table ##
#########################

t = get_halo_table(rs[0])

assert(len(t)>0)
assert(t.colnames[0]=='ID')

#############################
## test get_particle_table ##
#############################

p = get_particle_table(particles)

assert(len(p)>0)
assert(p.colnames[0]=='ID')

############################
## test special halo save ##
############################

h = Halo(0, rs, [i/10 for i in range(10)])
h.populate_full_particles(particles, False) 
h.save("halo0_sp.npy")

h2 = load_halo("halo0_sp.npy", "test_output/particles/")

assert(h==h2)

