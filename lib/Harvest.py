
import numpy as np
from astropy.table import Table 
from astropy.io import ascii

from Halo import *

# TODO: merge this usage with the halo infrastructure for harvest

def get_populated_halos(file_paths, rs_path):
    
    particles = get_particle_table(file_paths)
    halos = get_halo_table(rs_path)

    populated_halos = [] ; id_list = [] ; halo_index = 0
    for p in particles:
        halo_id = p['PHALO']
        while halo_index<len(halos) and halos['ID'][halo_index] != halo_id: halo_index += 1
        if halo_index==len(halos): break

        if len(populated_halos)==0 or populated_halos[-1].ID != halo_id: 
            populated_halos += [Halo()]
            populated_halos[-1].ID = halo_id
            
            populated_halos[-1].vx += [halos['VX'][halo_index]]
            populated_halos[-1].vy += [halos['VY'][halo_index]]
            populated_halos[-1].vz += [halos['VZ'][halo_index]]

        populated_halos[-1].particle_list += [Particle.Particle(p['X'],p['Y'],p['Z'],p['VX'],p['VY'],p['VZ'])]

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
            if len(data) < 10: continue
            x  += [data[0]] ; y  += [data[1]] ; z  += [data[2]]
            vx += [data[3]] ; vy += [data[4]] ; vz += [data[5]]
            IDs += [int(data[6])] ; pHalo += [int(data[9])]
        f.close()

    t = Table()
    t['ID'] = IDs ; t['PHALO'] = pHalo
    t['X']  = x  ; t['Y']  = y  ; t['Z']  = z
    t['VX'] = vx ; t['VY'] = vy ; t['VZ'] = vz

    try:
        return t.group_by('PHALO')
    except IndexError:
        print("WARNING: particle table is empty!")
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
        if line[:15] == "#Particle mass:":
            for n in line.split():
                try: m = float(n)
                except ValueError: continue

    return m

















