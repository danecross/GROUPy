
from Halo import load
from matplotlib import pyplot as plt

fig = plt.figure(figsize=(20.,20.))
ax = plt.axes(projection='3d')

h = load("halo_frommbp.npy", "test_output/MBP_save/")
t = 6

particles_x = [p.x[t] for p in h.particle_list]
particles_y = [p.y[t] for p in h.particle_list]
particles_z = [p.z[t] for p in h.particle_list]

ax.scatter3D(particles_x, particles_y, particles_z)


MBP = h.mbps[t]
ax.scatter([MBP.x[t]], [MBP.y[t]], [MBP.z[t]],color='r')

plt.savefig("maybe.png")


