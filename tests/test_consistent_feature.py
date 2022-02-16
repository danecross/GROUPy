
from Halo import *


files = ['test_data/really_consistent_%i.list'%i for i in range(10)]
times = [float(t) for t in open('test_data/output_snapshot_times.txt')]

#############################
## test halo intialization ## 
#############################

h = Halo(0, files, times, first_timestep=4, is_consistent=True, backtrack=False)

assert(type(h.PIDs)!=int)
assert(type(h.UPIDs)!=int)
assert(h.UPIDs[0] == -1 and h.UPIDs[5] == 0)
assert(h.PIDs[0] == -1 and h.PIDs[5] == 0)


###########################
## test get_halo_from_ID ##
###########################

h2 = h.get_halo_from_ID(1, 7, is_consistent=True)

assert(h.mass[5] == h2.mass[5])
assert(h != h2)
assert(h2.tag_idx == h.tag_idx)


#############################
## test get_nearest_parent ##
#############################

parent = h2.get_nearest_parent()

assert(parent == h)

#############################
## test get_largest_parent ##
#############################

ancestor = h2.get_largest_parent()

assert(ancestor == h)


