
import numpy as np
from astropy.table import Table 
from astropy.io import ascii

from Halo import *

def get_populated_halos(particle_paths, rs_path, timestep=0, num_timesteps=1):
    
    particles = get_particle_table(particle_paths)
    halos = get_halo_table(rs_path)

    populated_halos = [] ; id_list = [] ; halo_index = 0
    for p in particles:
        halo_id = p['PHALO']
        while halo_index<len(halos) and halos['ID'][halo_index] != halo_id: 
            halo_index += 1
        if halo_index==len(halos): break

        #check if we're putting the particle into the right halo
        if len(populated_halos)==0 or populated_halos[-1].ID != halo_id: 
            #make a new halo
            populated_halos += [Halo()]
            last_h = populated_halos[-1]
            last_h.ID = halo_id
            
            last_h.vx += [halos['VX'][halo_index]]
            last_h.vy += [halos['VY'][halo_index]]
            last_h.vz += [halos['VZ'][halo_index]]

        last_h.particle_list += [Particle.Particle(p['X'],p['Y'],p['Z'],\
						   p['VX'],p['VY'],p['VZ'],\
						   this_index=timestep, \
						   num_timesteps=num_timesteps)]

    return populated_halos


# gets the particle table from the rockstar outputs
# input:
#   file_paths: list of paths to the particle files
# output:
#   astropy Table of positions, velocities, and particle ID 
#       that is grouped by the parent halo id
def get_particle_table(file_paths):

    x  = [] ; y  = [] ; z  = []
    vx = [] ; vy = [] ; vz = []
    IDs = [] ; pHalo = []

    if type(file_paths) == str: file_paths = [file_paths]
    for path in file_paths:
        f = open(path,'r')
        for line in f:
            if line[0] == '#': continue
            data = [float(xi) for xi in line.split()]
            if len(data) < 9: continue
            x  += [data[0]] ; y  += [data[1]] ; z  += [data[2]]
            vx += [data[3]] ; vy += [data[4]] ; vz += [data[5]]
            IDs += [int(data[6])] ; pHalo += [int(data[-1])]
        f.close()

    t = Table()
    t['ID'] = IDs ; t['PHALO'] = pHalo
    t['X']  = x  ; t['Y']  = y  ; t['Z']  = z
    t['VX'] = vx ; t['VY'] = vy ; t['VZ'] = vz

    try:
        return t.group_by('PHALO')
    except IndexError:
        return t

# get and save the shortened particle table for one halo
# input: 
#   file_paths: list of paths to the particle file
#   halo_num: halo number to extract
#   save_path: path to save the file
# output: 
#   no return value, but saves file to save_path

def get_short_particle_table(file_paths, halo_num):

    x  = [] ; y  = [] ; z  = []
    vx = [] ; vy = [] ; vz = []
    IDs = [] ; pHalo = []

    if type(file_paths) == str: file_paths = [file_paths]
    for path in file_paths:
        f = open(path,'r')
        for line in f:
            if line[0] == '#': continue
            data = [float(xi) for xi in line.split()]
            #if len(data) < 9: continue
            if data[-1] != halo_num: continue
            x  += [data[0]] ; y  += [data[1]] ; z  += [data[2]]
            vx += [data[3]] ; vy += [data[4]] ; vz += [data[5]]
            IDs += [int(data[6])] ; pHalo += [int(data[-1])]
        f.close()

    t = Table()
    t['ID'] = IDs ; t['PHALO'] = pHalo
    t['X']  = x  ; t['Y']  = y  ; t['Z']  = z
    t['VX'] = vx ; t['VY'] = vy ; t['VZ'] = vz

    return t



# gets the table of halo information from rockstar outputs
# input:
#   rs_path: path to the rockstar out.*.list file
# output:
#   table of halo information ordered by ID
def get_halo_table(rs_path):

    hf = open(rs_path, 'r')
    header = hf.readline()[1:].split()
    hf.close()

    for i in range(len(header)): header[i]=header[i].upper()
    t = ascii.read(rs_path, format='no_header', names=header, delimiter=" ")
    t.sort('ID')
    return t

def get_particle_mass(rs_path):

    hf = open(rs_path, 'r')
    for line in hf:
        if "Particle mass" in line[:15]:
            for n in line.split():
                try: m = float(n)
                except ValueError: continue

    return m


# returns None if the line dne in the file
def get_time(rs_path):

    f = open(rs_path, 'r')
    a = None
    for line in f:
        if "a =" in line:
            for n in line.split():
                try: a = float(n)
                except ValueError: continue

    return a















