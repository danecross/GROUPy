
from Harvest import *
from Halo import *
from Particle import *
from Energies import *

out_path = "test_data/"
rs = [out_path+"out_%i.list"%i for i in range(10)]
particles = [out_path+"halos_snapshot_%i.hdf5.0.particles"%i for i in range(10)]

pdir = os.path.join("test_output", "particles")
num_timesteps = 10
pfile_base_stripped = "test_data/halos_snapshot_%i.hdf5"
all_pfiles_stripped = [pfile_base_stripped%i for i in range(num_timesteps)]
all_pfiles = [os.path.join(out_path, f) for f in os.listdir("test_data") if f[:10]=='halos_snap']
all_pfiles.sort()

f = open(out_path+"output_snapshot_times.txt", 'r')
times_list = [float(t) for t in f]

############################
## test harvest_particles ##
############################

for f in os.listdir(pdir):
    pf = os.path.join(pdir, f)
    os.remove(pf)

harvest_particles_by_timestep(all_pfiles[:2], 0, num_timesteps, pdir)

assert(len(os.listdir(pdir)) > 0)

for f in os.listdir(pdir):
    pf = os.path.join(pdir, f)
    os.remove(pf)

harvest_particles(all_pfiles, times_list, pdir)

assert(len(os.listdir(pdir)) > 0)

p = load_particle(os.path.join(pdir, "p_143.npy"))
assert((np.array(p.x)>-np.inf).all())


##############################
## test get particle masses ##
##############################

particle_mass = get_particle_mass(rs[0])
assert(particle_mass==4.61528e+06)

#########################
## test get_halo_table ##
#########################

t = get_halo_table(rs[0])

assert(len(t)>0)
assert(t.colnames[0]=='ID')

#############################
## test get_particle_table ##
#############################

p = get_particle_table(particles)

assert(len(p)>0)
assert(p.colnames[0]=='ID')

############################
## test special halo save ##
############################

h = Halo(0, rs, [i/10 for i in range(10)])
h.populate_full_particles(particles, False) 
h.save("halo0_sp.npy")

h2 = load_halo("halo0_sp.npy", "test_output/particles/")

assert(h==h2)

