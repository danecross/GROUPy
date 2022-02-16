

from Harvest import *

pfile = "test_data/halos_snapshot_2.hdf5.1.particles"

############################
## test get_particle_mass ##
############################

mass = get_particle_mass(pfile)

assert( mass == 4.61528e+06 )


###################
## test get_time ##
###################

a = get_time(pfile)

assert( a == 0.2 )

##############################
## test get_particle_bounds ##
##############################

bounds = get_particle_bounds(pfile)

assert(bounds[0] == (0.000000, 0.000000, 0.000000))
assert(bounds[1] == (3.105634, 2.702640, 0.970517))



