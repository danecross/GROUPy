
import os
from astropy.io import ascii
from astropy.table import Table
import numpy as np

from Particle import *
from Energies import *
from Harvest import *

import astropy

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

    def get_mb_particle(self, timestep):

        # input checking	
        if timestep > len(self.times_list): raise IndexError("desired timestep > number of available timesteps")
        if len(self.particle_list) == 0:
            raise Exception("NoParticlesAvailable","before finding the most bound particle, you must load the particle data using populate_full_particles or populate_particle_list")
        if len(self.mbps) == 0: self.mbps = [-1]*len(self.times_list)

        # if already calculated, return saved value
        if self.mbps[timestep] != -1: 
            pID = self.mbps[timestep]
            idx = np.where(np.array(self.particle_IDs) == pID)[0][0]
            return self.particle_list[idx]

        # else, calculate MBP
        halo_v = (self.vx[timestep], self.vy[timestep], self.vz[timestep])
        energies = [p.get_energy(timestep, self.particle_list, halo_v) for p in self.particle_list]
        energies = [e if e != np.nan else np.inf for e in energies]
        mbp_idx = np.where(energies == np.min(energies))[0][0] 
        MBP = self.particle_list[mbp_idx]

        self.mbps[timestep] = MBP.id
        return MBP

    def rank_particle_boundedness(self, timestep):

        if len(self.particle_list) == 0 : return
        halo_v = (self.vx[timestep], self.vy[timestep], self.vz[timestep])
        energies = [p.get_energy(timestep, self.particle_list, halo_v) for p in self.particle_list]
        ranked_particles = [p for _,p in sorted(zip(energies, self.particle_list))]
        
        if len(self.mbps) == 0: self.mbps = [-1]*len(self.times_list)
        self.mbps[timestep] = ranked_particles[0].id

        for p, i in zip(ranked_particles, range(len(ranked_particles))):
            p.bound_rank[timestep]=i

    def populate_particle_list(self, particle_file, timestep=None):

        if len(self.times_list)==0: raise Exception("EmptyHaloException", "Empty halo information; cannot make particle list")
        if self.ID == -1: raise Exception("InvalidHaloID","to read from particle table, must set the ID")
        if timestep is not None and timestep > len(self.ids): 
            raise Exception("InvalidTimestepIntertion", "attempting to insert particles beyond available timesteps")
        elif timestep is None:
            t0 = get_time(particle_file) 
            timestep =  np.where(np.abs(t0-np.array(self.times_list))<2e-5)[0][0]

        p_table = get_short_particle_table(particle_file, self.ids[timestep]) 
        if len(p_table)==0: return

        i = 0
        for tp in p_table:
            if tp['PHALO'] != self.ids[timestep]: continue
            idx = np.where(self.particle_IDs==tp['ID'])[0]
            if len(idx) > 0: # particle is already in the list
                p = self.particle_list[idx[0]] 
            else: # a new particle
                p = Particle(num_timesteps=len(self.times_list), id=tp['ID'],\
				this_index=timestep, mass=get_particle_mass(particle_file))
                self.particle_IDs += [tp['ID']]
                self.particle_list += [p]; i+=1

            p.insert_at_timestep(tp['X'],tp['Y'],tp['Z'],tp['VX'],tp['VY'],tp['VZ'], timestep, parent_halo=self.ID)


    def populate_full_particles(self, particles_list, find_all_files=True):

        if find_all_files: files_list = [self._get_p_names_list(p_name) for p_name in particles_list]
        else: files_list = [[f,] for f in particles_list]

        for files in files_list:
            t0 = get_time(files[0])
            t_idxs = np.where(np.abs(t0-np.array(self.times_list))<2e-5)[0]
            if len(t_idxs) == 0:
                raise Exception("InvalidTimesList","file is at timestep a=%f, no match in times_list. file name: %s"%(t0,files[0]))
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

    def load_particle(self, pid, pdir=None):
    
        p_file = os.path.join(pdir, "p_%i.npy"%pid)
        if not os.path.exists(p_file): raise FileNotFoundError("file %s does not exist"%p_file)
        
        p = load_particle(p_file) ; is_parent = False
        for ID, pid in zip(self.ids, p.parent_ids): 
            if int(ID) == int(pid): is_parent=True
        if not is_parent: raise Exception("InvalidParticlePlacing", "particle %i is never apart of halo %i"%(p.id, self.ID))

        self.particle_list += [p]

    def get_desc_list(self, first_id, files, times, file_start=0, verbose=False): 

        self.ids = np.array([-1]*(file_start))
        self.mass = np.array([-1]*(file_start))
        self.radius = np.array([-1]*(file_start))
        self.t = np.array([-1]*(file_start))

        self.x = np.array([-1]*(file_start))
        self.y = np.array([-1]*(file_start))
        self.z = np.array([-1]*(file_start))
        
        self.vx = np.array([-1]*(file_start)) 
        self.vy = np.array([-1]*(file_start)) 
        self.vz = np.array([-1]*(file_start))

        if self.ID == -1: self.ID = first_id
        if verbose: print("\ttracking halo %i"%self.ID)
        i = file_start ; last_id = first_id
        
        hf = open(files[0], 'r')
        header = hf.readline()[1:].split()
        hf.close()

        for f,ti in zip(files[file_start:], times[file_start:]):

            try:
                t = ascii.read(f, format='no_header', names=header, delimiter=" ")
            except astropy.io.ascii.core.InconsistentTableError:
                self.add_timestep(i, -1, -1, -1, -1, -1, -1, -1, -1, -1, ti)
                continue

            # get row
            ri = np.where(t['ID']==last_id)[0]

            if len(ri) == 0: 
                self.add_timestep(i, -1, -1, -1, -1, -1, -1, -1, -1, -1, ti)
                continue

            self.add_timestep(i, last_id, t['Mvir'][ri], t['Rvir'][ri], \
                                 t['X'][ri], t['Y'][ri], t['Z'][ri], \
                                 t['VX'][ri], t['VY'][ri], t['VZ'][ri], ti)
            
            i += 1 ; last_id = t['DescID'][ri]

            if verbose and (i-file_start)%50==0: 
                print("\t\ttracked halo %i through %i screenshots "%(self.ID,i))

            if last_id == -1: # the halo has no descendant (it disappeared) 
                while i < len(files):
                    self.add_timestep(i, -1, -1, -1, -1, -1, -1, -1, -1, -1, ti)
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

    def get_distances_from(self, x0=None, y0=None, z0=None, other=None):

        if other is not None: self._get_distances_from(other)

        x = np.array([xi-x0 for xi in self.x if xi != -1])
        y = np.array([yi-y0 for yi in self.y if yi != -1])
        z = np.array([zi-z0 for zi in self.z if zi != -1])

        distances = np.sqrt(x**2+y**2+z**2)

        return distances

    def _get_distances_from(self, other):

        x = np.array([xi-x0 for xi, x0 in zip(self.x, other.x) if xi != -1 and x0 != -1])
        y = np.array([yi-y0 for yi, y0 in zip(self.y, other.y) if yi != -1 and y0 != -1])
        z = np.array([zi-z0 for zi, z0 in zip(self.z, other.z) if zi != -1 and z0 != -1])

        distances = np.sqrt(x**2+y**2+z**2)

        return distances

    # particle dir must exist before running save
    def save(self, name="halo.npy"):

        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) \
					     and not attr.startswith("__")\
					     and not attr=="particle_list"]

        h = {attr:getattr(self,attr) for attr in members}

        np.save(name, h, allow_pickle=True)
 
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
#   halo object loaded from given file
def load_halo(name, particle_directory=None):

    h = Halo()
    h_dict = np.load(name, allow_pickle=True)

    for attr in h_dict.item().keys():
        if attr == "particle_list": continue
        setattr(h,attr, h_dict.item()[attr])

    if particle_directory is None: return h
    if len(h.particle_IDs) == 0: 
        raise Exception("NoParticlesFound", "halo's list of member particle IDs is empty")
    
    for pid in h.particle_IDs:
        h.load_particle(pid, pdir=particle_directory)

    return h












