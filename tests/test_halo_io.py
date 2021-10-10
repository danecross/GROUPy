
import os 
import numpy as np
from astropy.io import ascii

from Halo import Halo, load

#################################
## test empty halo constructor ##
#################################

empty_halo = Halo()

assert(hasattr(empty_halo, "ids") and len(empty_halo.ids)==0)

#######################################################
## test case for non-empty constructor: no particles ##
#######################################################

first_idx = 4
base = "test_data/out_%i.list"
file = base%first_idx
files = [base%i for i in range(100) if os.path.exists(base%i)]

times_file = open("test_data/output_snapshot_times.txt", 'r')
times = [float(line) for line in times_file]

###############################
## initialize non-empty halo ##
###############################

f = open(file, 'r')
header = f.readline()[1:].split()
f.close()
t = ascii.read(file, format='no_header', names=header, delimiter=" ")

line = t[0]
halo = Halo(line['ID'], files, times, first_index=first_idx)

assert(hasattr(halo, "ids")); assert(hasattr(halo, "mass")) ; assert(hasattr(halo, "radius"))
assert(hasattr(halo, "x")); assert(hasattr(halo, "y")) ; assert(hasattr(halo, "z"))
assert(hasattr(halo, "vx")); assert(hasattr(halo, "vy")) ; assert(hasattr(halo, "vz"))

# test that backtrack works
assert(halo.x[first_idx-1]!=-1)


#####################################################
## initialize new halo with fully loaded particles ##
#####################################################

particle_files = ["test_data/halos_snapshot_%i.hdf5"%i for i in range(10)]

f = open(file, 'r')
header = f.readline()[1:].split()
f.close()
t = ascii.read(file, format='no_header', names=header, delimiter=" ")

line = t[0]
halo = Halo(line['ID'], files, times, particle_files, first_index=first_idx)

assert(hasattr(halo, "particle_list")) 
assert(len(halo.particle_list) > 0)

p = halo.particle_list[0]
assert(hasattr(p, "x")) 
for x in p.x:
    assert(x > -np.inf)

##############################
## test save/load procedure ##
##############################

test_dir = "test_output/"
if not os.path.exists(test_dir): os.mkdir(test_dir)

halo.save(name="halo%i.npy"%halo.ID, particle_dir=test_dir+"halo_%i"%halo.ID)
assert(os.path.exists(test_dir+"halo_0/halo%i.npy"%halo.ID))

halo2 = load(test_dir+"halo_0/halo%i.npy"%0)

assert(halo==halo2)

