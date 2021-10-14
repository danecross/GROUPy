
from Halo import *

import os
from matplotlib import pyplot as plt


### an example on how to create a halo from the out_*.list files

## using the data from the tests..
rs_base = "../../tests/test_data/out_%i.list"
halo_data = [rs_base%i for i in range(10)]

time_file = "../../tests/test_data/output_snapshot_times.txt"
tf = open(time_file, 'r')
times = [float(l) for l in tf]


## create a halo (no particle data input)

h = Halo(0, halo_data, times)

## plot some data. 
### note: these plots will be rather boring with the data from the test
###       directory, as the mass and radii don't change as a function of time.
###       for more fun, add data from actual simulations

output_dir = "plotting_output/"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


### mass as a function of time

plt.plot(times, h.mass)
plt.savefig(os.path.join(output_dir,"mass.png"))
plt.cla()

### radius as a function of time

plt.plot(times, h.radius)
plt.savefig(os.path.join(output_dir,"radius.png"))


