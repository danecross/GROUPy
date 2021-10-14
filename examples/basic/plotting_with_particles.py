
from Halo import *

import os
from matplotlib import pyplot as plt


### an example on how to create a halo from the out_*.list files

## using the data from the tests..

pd_base = "../../tests/test_data/halos_snapshot_%i.hdf5"
particle_data = [pd_base%i for i in range(10)]

rs_base = "../../tests/test_data/out_%i.list"
halo_data = [rs_base%i for i in range(10)]

time_file = "../../tests/test_data/output_snapshot_times.txt"
tf = open(time_file, 'r')
times = [float(l) for l in tf]


## create a halo (no particle data input)

h = Halo(0, halo_data, times, particle_data)

## plot some data. 
### note: these plots will be rather chaotic with the data from the test
###       directory, as the data is randomly generated.
###       for more fun, add data from actual simulations

output_dir = "plotting_output/"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

### energy as a function of time

p = h.get_mb_particle(3)
[p.get_energy(i, h.particle_list) for i in range(len(times))]

plt.plot(times, p.energy)
plt.savefig(os.path.join(output_dir,"p_energy.png"))
plt.cla()

### boundedness as a function of time

[h.rank_particle_boundedness(i) for i in range(len(times))]

plt.plot(times, p.bound_rank)
plt.savefig(os.path.join(output_dir,"bound_rank.png"))


