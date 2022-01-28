
from Particle import *
from Halo import *
from Harvest import *


num_timesteps = 10
pfile_base_stripped = "test_data/halos_snapshot_%i.hdf5"
all_pfiles_stripped = [pfile_base_stripped%i for i in range(num_timesteps)]
all_pfiles = [f for f in os.listdir("test_data") if f[:10]=='halos_snap']

##################################
## test basic particle creation ##
##################################

# Empty Particle Init
# as long as this doesn't throw and error, success
p = Particle()

# Non-empty Particle Init

p = Particle(1, 1, 1, 1, 1, 1)
assert(len(p.x)==1)

p = Particle(1, 1, 1, 1, 1, 1, this_index=1, num_timesteps=2)
assert(len(p.x)==2)

p = Particle(1, 1, 1, 1, 1, 1, mass=1)
assert(p.mass==1)

#############################
## test insert_at_timestep ##
#############################

p = Particle(1, 1, 1, 1, 1, 1, this_index=1, num_timesteps=5)
p.insert_at_timestep(2, 2, 2, 2, 2, 2, 3, 123)
assert(p.x[3]==2)


###############
## test save ##
###############

pdir = "test_output/basic_pout/"
pid = 20
if not os.path.exists(pdir): os.mkdir(pdir)
if os.path.exists(pdir+"p_%i.npy"%pid): os.remove(pdir+"p_%i.npy"%pid)

p = Particle(1, 1, 1, 1, 1, 1, id=pid)
p.save(pdir)

assert(os.path.exists(pdir+"p_%i.npy"%pid))


########################
## test load_particle ##
########################

pl = load_particle(pdir+"p_%i.npy"%pid)

assert(pl==p)
assert(len(pl.x) == len(p.x))














