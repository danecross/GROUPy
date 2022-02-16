
import os
import shutil

from Particle import *
from Halo import * 
from Harvest import *

###################
## get test data ##
###################

pfiles = ["test_data/halos_snapshot_%i.hdf5"%i for i in range(10)]
hfiles = ["test_data/out_%i.list"%i for i in range(10) if os.path.exists("test_data/out_%i.list"%i)]

f = open('test_data/output_snapshot_times.txt','r')
times = [float(line) for line in f]
f.close()

###############
## make halo ## 
###############

h = Halo(0, hfiles, times, pfiles, first_timestep=4, backtrack=False)


####################################
## find most bound particle ranks ##
####################################

for i in range(len(times)):
    h.rank_particle_boundedness(i)

p = h.particle_list[0]
assert(all(rank >= 0 for rank in p.bound_rank))

##############################
## find most bound particle ##
##############################

MBP = h.get_mb_particle(6)

assert(all(rank >= 0 for rank in MBP.bound_rank))




