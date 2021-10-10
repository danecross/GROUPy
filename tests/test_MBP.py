
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

savedir = "test_output/MBP_save/"
if os.path.exists(savedir):
    shutil.rmtree(savedir)
    os.mkdir(savedir)

############################################
## make halo and find most bound particle ##
############################################

h = Halo(0, hfiles, times, pfiles, first_timestep=4, backtrack=False)
h.get_mb_particle(6)

h.save("halo_frommbp.npy", "test_output/MBP_save/")








