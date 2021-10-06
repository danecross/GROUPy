
import numpy as np
from matplotlib import pyplot as plt
from Halo import *


## using the data from the tests..

pd_base = "../../tests/test_data/halos_snapshot_%i.hdf5.0.particles"
particle_data = [pd_base%i for i in range(10)]

rs_base = "../../tests/test_data/out_%i.list"
halo_data = [rs_base%i for i in range(10)]

time_file = "../../tests/test_data/output_snapshot_times.txt"
tf = open(time_file, 'r')
times = [float(l) for l in tf]

## ..find the most bound particle for timestep 4

t = 4
h = Halo(0, halo_data, times, particle_data, file_start=4)

MBP = h.get_mb_particle(4)

## then graph some stuff!

# first, let's graph what the data looks like in 3D:

fig = plt.figure(figsize=(20.,20.))
ax = plt.axes(projection='3d')

particles_x = [p.x[t] for p in h.particle_list]
particles_y = [p.y[t] for p in h.particle_list]
particles_z = [p.z[t] for p in h.particle_list]

vels = [np.sqrt(p.x[t]**2+p.y[t]**2+p.z[t]**2) for p in h.particle_list]

ax.scatter3D(particles_x, particles_y, particles_z, s=200, c=vels)

MBP = h.mbps[t]
ax.scatter([MBP.x[t]], [MBP.y[t]], [MBP.z[t]],facecolor=(0,0,0,0),s=1000, edgecolor='red')

plt.savefig("3d-render.png")




