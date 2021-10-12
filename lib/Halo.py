
import os
from astropy.io import ascii
from astropy.table import Table
import numpy as np

import Particle 
from Energies import *
from Harvest import *

class Halo(object):

    def __init__(self, *args, **kwargs): 

        first_idx = kwargs.get('first_timestep', 0)
        backtrack = kwargs.get('backtrack', True)
        verbose = kwargs.get('verbose', False) 

        self.halo_files = [] 

        self.mbps = []
        self.particle_list = [] ; self.particle_IDs = []
        self.ID = -1

        if len(args) == 0:
            self.ids = [] ; self.tag_idx = -1
            self.mass = []
            self.radius = []
            self.x = [] ; self.y = [] ; self.z = [] 
            self.vx = [] ; self.vy = [] ; self.vz = [] 
        elif len(args) >= 3:

            tag_id = args[0] 
            self.halo_files = args[1]
            self.times_list = args[2]
            
            self.tag_idx = first_idx
            self.get_desc_list(tag_id, self.halo_files, self.times_list, file_start=first_idx, verbose=verbose)
            if backtrack:
                self.backtrack(tag_id, self.halo_files, self.times_list, file_start=first_idx)
            if len(args)==4: 
                particle_files = args[3]
                self.populate_full_particles(particle_files)

        else:
            raise TypeError("Two options for initializing a Halo object: \n" +\
                            "\t 1. default constructor (no arguments) \n " +\
                            "\t 2. full constructor w/ 3 required arguments and 3 optional keyword arguments: \n" +\
                            "\t\t ID of the first halo appearance in the rockstar files (int),\n " +\
                            "\t\t time-ordered list of rockstar files (list)\n " +\
                            "\t\t list of times for each snapshot (t = 1/(1+z)) (list)\n " +\
                            "\t\t index of file where the halo first appears (int, optional, default=0)\n"+\
                            "\t\t option to backtrack halo information from the tagging redshift (bool, optional, default=True)\n"+\
                            "\t\t flag for verbose output (bool, optional, default=False)\n"+\
                            "\t 3. full constructor with particle population\n"+\
                            "\t\t ID of the first halo appearance in the rockstar files (int),\n " +\
                            "\t\t time-ordered list of rockstar files (list)\n " +\
                            "\t\t list of times for each snapshot (t = 1/(1+z)) (list)\n " +\
                            "\t\t array of particle path names ordered by time\n"+\
                            "\t\t index of file where the halo first appears (int, optional, default=0)\n"+\
                            "\t\t option to backtrack halo information from the tagging redshift (bool, optional, default=True)\n"+\
                            "\t\t flag for verbose output (bool, optional, default=False)\n"+\
                            "len of args: %i, len of kwargs: %i"%(len(args), len(kwargs)))
            exit()

    def get_mb_particle(self, timestep):

        # input checking	
        if timestep > len(self.times_list): raise IndexError("desired timestep > number of available timesteps")
        if len(self.particle_list) == 0 or self.particle_list[0].x[timestep] == -np.inf : 
            raise Exception("NoParticlesAvailable","before finding the most bound particle, you must load the particle data using populate_full_particles or populate_particle_list")
        if len(self.mbps) == 0: self.mbps = [-1]*len(self.times_list)

        # if already calculated, return saved value
        if self.mbps[timestep] != -1: return self.mbps[timestep]

        # else, calculate MBP
        halo_v = (self.vx[timestep], self.vy[timestep], self.vz[timestep])
        energies = [p.get_energy(timestep, self.particle_list, halo_v) for p in self.particle_list]
        mbp_idx = np.where((energies == np.min(energies)))[0][0] 
        MBP = self.particle_list[mbp_idx]

        self.mbps[timestep] = MBP
        return MBP

    def rank_particle_boundedness(self, timestep):

        halo_v = (self.vx[timestep], self.vy[timestep], self.vz[timestep])
        energies = [p.get_energy(timestep, self.particle_list, halo_v) for p in self.particle_list]
        ranked_particles = [p for _,p in sorted(zip(energies, self.particle_list))]
        for p, i in zip(ranked_particles, range(len(ranked_particles))):
            p.bound_rank[timestep]=i

    def populate_particle_list(self, particle_file, timestep):
	
        if len(self.times_list)==0: raise Exception("EmptyHaloException", "Empty halo information; cannot make particle list")
        if self.ID == -1: raise Exception("InvalidHaloID","to read from particle table, must set the ID")

        p_table = get_particle_table(particle_file)
        mask = (p_table['PHALO']==self.ID)
        for tp in p_table[mask]:
            idx = np.where(self.particle_IDs==tp['ID'])[0]
            if len(idx) > 0: # particle is already in the list
                p = self.particle_list[idx[0]]
            else: # a new particle
                p = Particle.Particle(num_timesteps=len(self.times_list),\
				this_index=timestep, mass=get_particle_mass(particle_file))
                self.particle_IDs += [tp['ID']]
                self.particle_list += [p]

            p.insert_at_timestep(tp['X'],tp['Y'],tp['Z'],tp['VX'],tp['VY'],tp['VZ'], timestep)
				

    def populate_full_particles(self, particles_list):

        files_list = [self._get_p_names_list(p_name) for p_name in particles_list]

        for files in files_list:
            t0 = get_time(files[0])
            t_idxs = np.where(t0==np.array(self.times_list))[0]
            if len(t_idxs) == 0:
                raise Exception("InvalidTimesList","file is at timestep a=%.2f, no match in times_list"%t0)
            t_idx = t_idxs[0]
            for f in files:
                # check all files have same time 
                t = get_time(f)
                if t != t0: 
                    raise Exception("InconsistentFiles","file %s has different a-value (time) than the first file")
                
                self.populate_particle_list(f, t_idx)


    def _get_p_names_list(self, pfile_name):
        
        base = pfile_name+".%i.particles"
        res = [base%i for i in range(0,1000) if os.path.exists(base%i)]
    
        if len(res) == 0: 
            raise Exception("InvalidFileBase","could not find any particle files %s"%(base))

        return res


    def get_desc_list(self, first_id, files, times, file_start=0, verbose=False): 

        self.ids = np.array([-1]*(file_start))
        self.mass = np.array([-1]*(file_start))
        self.radius = np.array([-1]*(file_start))
        self.t = np.array([-1]*(file_start))
        self.x = np.array([-1]*(file_start))  ; self.y = np.array([-1]*(file_start))  ; self.z = np.array([-1]*(file_start))
        self.vx = np.array([-1]*(file_start)) ; self.vy = np.array([-1]*(file_start)) ; self.vz = np.array([-1]*(file_start))

        if self.ID == -1: self.ID = first_id
        if verbose: print("\ttracking halo %i"%self.id)
        i = file_start ; last_id = first_id
        
        hf = open(files[0], 'r')
        header = hf.readline()[1:].split()
        hf.close()

        for f,ti in zip(files[file_start:], times[file_start:]):

            t = ascii.read(f, format='no_header', names=header, delimiter=" ")

            # get row
            ri = np.where(t['ID']==last_id)[0]

            self.add_timestep(i, last_id, t['Mvir'][ri], t['Rvir'][ri], \
                                 t['X'][ri], t['Y'][ri], t['Z'][ri], \
                                 t['VX'][ri], t['VY'][ri], t['VZ'][ri], ti)
            
            i += 1 ; last_id = t['DescID'][ri]

            if verbose and (i-file_start)%50==0: print("\t\ttracked halo %i through %i screenshots "%(self.id,i))

            if last_id == -1: # the halo has no descendant (it disappeared) 
                while i < len(files):
                    self.add_timestep(i, -1, -1, -1, -1, -1, -1, -1, -1, -1, t)
                    i += 1
                break

    def backtrack(self, first_id, files, times, file_start=0, verbose=False):
        
        this_id = first_id ; i = file_start-1
        for f,ti in zip(reversed(files[:file_start]), reversed(times[:file_start])):
            
            hf = open(f, 'r')
            header = hf.readline()[1:].split()
            hf.close()
            t = ascii.read(f, format='no_header', names=header, delimiter=" ")

            try:
                ri = np.where(t['DescID']==this_id)[0][0]
            except IndexError:
                return

            self.add_timestep(i, this_id, t['Mvir'][ri], t['Rvir'][ri], \
                                 t['X'][ri], t['Y'][ri], t['Z'][ri], \
                                 t['VX'][ri], t['VY'][ri], t['VZ'][ri], ti)

            i -= 1 ; this_id = t['ID'][ri]

    def add_timestep(self, index, ID, mass, radius, x, y, z, vx, vy, vz, t):

        if index == len(self.ids):
            self.ids	 = np.append(self.ids	 , ID)
            self.mass    = np.append(self.mass   , mass)
            self.radius  = np.append(self.radius , radius)
            self.x       = np.append(self.x      , x)
            self.y       = np.append(self.y      , y)
            self.z       = np.append(self.z      , z)
            self.vx      = np.append(self.vx     , vx)
            self.vy      = np.append(self.vy     , vy)
            self.vz      = np.append(self.vz     , vz)
            self.t       = np.append(self.t      , t)
        elif index < len(self.ids) and self.ids[index]==-1:
            self.ids[index]	= ID
            self.mass[index]    = mass
            self.radius[index]  = radius
            self.x[index]       = x
            self.y[index]       = y 
            self.z[index]       = z 
            self.vx[index]      = vx 
            self.vy[index]      = vy
            self.vz[index]      = vz
            self.t[index]       = t
        elif index < len(self.ids) and self.ids[index]!=-1:
            raise ValueError("attempting to change read values from the files")
            exit()
        else:
            raise ValueError("attempting to add a timestep which has no previous timestep")

    def get_distances_from(self, x0, y0, z0):

        x = np.array([xi-x0 for xi in self.x if xi != -1])
        y = np.array([yi-y0 for yi in self.y if yi != -1])
        z = np.array([zi-z0 for zi in self.z if zi != -1])

        distances = np.sqrt(x**2+y**2+z**2)

        return distances

    def get_distances_from(self, other):

        x = np.array([xi-x0 for xi, x0 in zip(self.x, other.x) if xi != -1 and x0 != -1])
        y = np.array([yi-y0 for yi, y0 in zip(self.y, other.y) if yi != -1 and y0 != -1])
        z = np.array([zi-z0 for zi, z0 in zip(self.z, other.z) if zi != -1 and z0 != -1])

        distances = np.sqrt(x**2+y**2+z**2)

        return distances

    def save(self, name="halo.npy", particle_dir="halo/"):
        
        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) \
					     and not attr.startswith("__")\
					     and not attr=="particle_list"]

        h = {attr:getattr(self,attr) for attr in members}
	
        if len(self.particle_list) == 0: 
            np.save(name, h, allow_pickle=True)
            return 

        if not os.path.exists(particle_dir):
            os.mkdir(particle_dir)

        for p,i in zip(self.particle_list, range(len(self.particle_list))):
            p.save(os.path.join(particle_dir,"p_%05i.npy"%i))

        np.save(os.path.join(particle_dir,name), h, allow_pickle=True)

    def __str__(self):
        num_ss = len([x for x in self.ids if x != -1])
        return "Halo %i evolution with %i timesteps"%(self.ID, num_ss)

    def __eq__(self, other):
        
        if not isinstance(other, self.__class__): return False
        return (self.ID == other.ID and self.ID != -1 and other.ID != -1)


# loads a saved halo from directory
# input:
#   name: name of the file that contains the saved halo
# output; 
#   halo object that was saved
def load(name, particle_directory=None):

    h = Halo()
    
    if particle_directory is not None: fname = os.path.join(particle_directory, name)
    else: fname = name
    
    h_dict = np.load(fname, allow_pickle=True)

    for attr in h_dict.item().keys():
        if attr == "particle_list": continue
        setattr(h,attr, h_dict.item()[attr])

    if particle_directory is None: return h
    
    for fname in os.listdir(particle_directory):
        p = Particle.load(os.path.join(particle_directory,fname))
        h.particle_list += [p]

    return h












