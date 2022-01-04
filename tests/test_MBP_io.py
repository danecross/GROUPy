
from Harvest import *
from Halo import *
from Energies import *

out_path = "test_data/"
rs = out_path+"out_4.list"
particles = [out_path+"halos_snapshot_4.hdf5.0.particles"]


##############################
## test get particle masses ##
##############################

particle_mass = get_particle_mass(rs)
assert(particle_mass==4.61528e+06)

#########################
## test get_halo_table ##
#########################

t = get_halo_table(rs)

assert(len(t)>0)
assert(t.colnames[0]=='ID')

#############################
## test get_particle_table ##
#############################

p = get_particle_table(particles)

assert(len(p)>0)
assert(p.colnames[0]=='ID')

##############################
## test get_populated_halos ##
##############################

halos = get_populated_halos(particles, rs)

assert(type(halos)==list)
assert(type(halos[0])==Halo)

############################
## test special halo save ##
############################

h = halos[0]
h.save("halo0_sp.npy", "test_output/io")

h2 = load_halo("halo0_sp.npy", "test_output/io")

assert(h==h2)

