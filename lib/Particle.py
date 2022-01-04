
import numpy as np

from Energies import *

class Particle(object):

    def __init__(self, *args, **kwargs):

        this_idx = kwargs.get('this_index', 0)
        num_timesteps = kwargs.get('num_timesteps', 1)
        self.mass = kwargs.get('mass',-np.inf)
        self.id = kwargs.get('id',-1)

        if num_timesteps < 1: raise Exception("negative timesteps","number of timesteps must be >=1")
        elif this_idx > num_timesteps:
            raise Exception("this_index > num_timesteps", "insertion index is beyond the number of timesteps")

        self.x  = [-np.inf]*num_timesteps ; self.vx = [-np.inf]*num_timesteps
        self.y  = [-np.inf]*num_timesteps ; self.vy = [-np.inf]*num_timesteps
        self.z  = [-np.inf]*num_timesteps ; self.vz = [-np.inf]*num_timesteps
        self.energy = [-np.inf]*num_timesteps
        self.bound_rank = [-1]*num_timesteps

        self.member_of = [-1]*num_timesteps

        if len(args) == 6:
            i = this_idx
            self.x[i]  = args[0] ; self.y[i]  = args[1] ; self.z[i]  = args[2]
            self.vx[i] = args[3] ; self.vy[i] = args[4] ; self.vz[i] = args[5]
       

    def insert_at_timestep(self, x, y, z, vx, vy, vz, timestep, parent_halo=-1):

        if timestep > len(self.x): raise IndexError("timestep > number of timesteps")
        self.x[timestep]  = x  ; self.y[timestep]  = y  ; self.z[timestep]  = z
        self.vx[timestep] = vx ; self.vy[timestep] = vy ; self.vz[timestep] = vz

        self.member_of[timestep] = parent_halo


    def get_energy(self, timestep, particle_field, halo_v=(0,0,0)):

        if self.energy[timestep] > -np.inf: return self.energy[timestep]
        if self.mass == -np.inf: raise ValueError("must set particle mass before calculating energies") 
        K = KE(self, timestep, halo_v)
        P = PEg(self, particle_field, timestep) 

        self.energy[timestep] = K+P
        return K+P

    def save(self, name="particle.npy"):

        members = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]

        p = {attr:getattr(self,attr) for attr in members}
        np.save(name, p, allow_pickle=True)


    def __str__(self):
        return "particle %i"%self.id

    def __eq__(self, other):

        if not isinstance(other, self.__class__): return False
        return (self.id == other.id and self.id != -1 and other.id != -1)

def load_particle(name):

    p_dict = np.load(name, allow_pickle=True)
    p = Particle(num_timesteps=len(p_dict.item().get('x')))

    for attr in p_dict.item().keys():

        setattr(p,attr, p_dict.item()[attr])

    return p






