
import numpy as np

# returns potential energy in units of Msun*Mpc^2/s^2
def PEg(p, particle_field, timestep):

    G = 4.5155e-48                  # Mpc^3/Msun/s^2
    particle_mass = p.mass
    if particle_mass == -np.inf: raise ValueError("must set particle mass before calculating energies")

    PE = 0 ; i = timestep
    for op in particle_field:
        r = np.sqrt((p.x[i]-op.x[i])**2 + (p.y[i]-op.y[i])**2 + (p.z[i]-op.z[i])**2)
        if r == 0: continue # skip self-field
        PE += -G*particle_mass**2/r

    return PE

def KE(p, timestep, halo_v=(0,0,0)):
    
    particle_mass = p.mass
    if particle_mass == -np.inf: raise ValueError("must set particle mass before calculating energies")
    
    vx = (p.vx[timestep]-halo_v[0])/3.086e19             # convert from km/s to Mpc/s
    vy = (p.vy[timestep]-halo_v[1])/3.086e19
    vz = (p.vz[timestep]-halo_v[2])/3.086e19
    v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
    return 0.5*particle_mass*v_mag**2



