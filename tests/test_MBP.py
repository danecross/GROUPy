
import os
import shutil

from Particle import *
from Halo import * 
from Harvest import *

# get test data
pfiles = ["test_data/halos_snapshot_%i.hdf5.0.particles"%i for i in range(10)]
hfiles = ["test_data/out_%i.list"%i for i in range(10) if os.path.exists("test_data/out_%i.list"%i)]

f = open('test_data/output_snapshot_times.txt','r')
times = [float(line) for line in f]
f.close()

savedir = "test_output/MBP_save/"
if os.path.exists(savedir):
    shutil.rmtree(savedir)
    os.mkdir(savedir)

h = Halo(0, hfiles, times, first_timestep=4, backtrack=False)
h.particle_files = pfiles 
h.get_mb_particle(6)

h.save("halo_frommbp.npy", "test_output/MBP_save/")








