

from numpy.random import normal
from Particle import *

preamble = "#ID DescID Mvir Rvir X Y Z VX VY VZ \n\
#a = %.6f \n\
#Om = 0.279000; Ol = 0.721000; h = 0.700000 \n\
#FOF linking length: 0.280000 \n\
#Unbound Threshold: 0.500000; FOF Refinement Threshold: 0.700000 \n\
#Particle mass: 4.61528e+06 Msun/h \n\
#Box size: 10.000000 Mpc/h \n\
#Force resolution assumed: 0.005 Mpc/h \n\
#Units: Masses in Msun / h \n\
#Units: Positions in Mpc / h (comoving) \n\
#Units: Velocities in km / s (physical, peculiar) \n\
#Units: Halo Distances, Lengths, and Radii in kpc / h (comoving) \n\
#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical) \n\
#Units: Spins are dimensionless \n\
#Np is an internal debugging quantity. \n\
#Rockstar Version: 0.99.9-RC3+\n"

p_preamble = "#Particle table:\n\
#x y z vx vy vz particle_id assigned_internal_haloid internal_haloid external_haloid\n\
#Notes: As not all halos are printed, some halos may not have external halo ids.  (Hence the need to print internal halo ids).  Each particle is assigned to a unique halo; however, some properties (such as halo bound mass) are calculated including all substructure.  As such, particles belonging to subhalos are included in outputs; to exclude substructure, verify that the internal halo id is the same as the assigned internal halo id.\n\
#a = %.6f\n\
#Bounds: (0.000000, 0.000000, 0.000000) - (3.105634, 2.702640, 0.970517)\n\
#Om = 0.279000; Ol = 0.721000; h = 0.700000\n\
#FOF linking length: 0.280000\n\
#Unbound Threshold: 0.500000; FOF Refinement Threshold: 0.700000\n\
#Particle mass: 4.61528e+06 Msun/h\n\
#Box size: 10.000000 Mpc/h\n\
#Total particles processed: 349530\n\
#Force resolution assumed: 0.005 Mpc/h\n\
#Units: Masses in Msun / h\n\
#Units: Positions in Mpc / h (comoving)\n\
#Units: Velocities in km / s (physical, peculiar)\n\
#Units: Halo Distances, Lengths, and Radii in kpc / h (comoving)\n\
#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n\
#Units: Spins are dimensionless\n\
#Units: Total energy in (Msun/h)*(km/s)^2 (physical)\n\
#Note: idx, i_so, and i_ph are internal debugging quantities\n\
#Np is an internal debugging quantity.\n\
#Rockstar Version: 0.99.9-RC3+\n"


# Generate data in 10 timesteps with 20 particles of mass 1 Msun
num_timesteps = 10 ; num_pfiles = 2
num_p = 200

times_list = [i/num_timesteps for i in range(num_timesteps)]

h_ID = [0]*num_timesteps ; h_DescID = [0]*num_timesteps
h_X = normal(loc=0, scale=1, size=num_timesteps)
h_Y = normal(loc=0, scale=1, size=num_timesteps)
h_Z = normal(loc=0, scale=1, size=num_timesteps)
h_VX = normal(loc=0, scale=1, size=num_timesteps)
h_VY = normal(loc=0, scale=1, size=num_timesteps)
h_VZ = normal(loc=0, scale=1, size=num_timesteps)

Mvir = [num_p]*num_timesteps ; Rvir = [1]*num_timesteps

###################
# write halo data #
###################

base_rs = "out_%i.list"
rs_files = [base_rs%i for i in range(num_timesteps)]

for fn, i in zip(rs_files,range(num_timesteps)):
    
    f = open(fn,'w')

    l = preamble%(times_list[i])
    f.write(l)

    line = [h_ID[i], h_DescID[i], Mvir[i], Rvir[i], \
	    h_X[i], h_Y[i], h_Z[i], h_VX[i], h_VY[i], h_VZ[i]]
    line = [str(l) for l in line]

    f.write(" ".join(line))

    f.close()

#######################
# write particle data #
#######################

particles = [Particle(num_timesteps=num_timesteps, mass=1) for _ in range(num_p)]

ids = [i for i in range(num_p)]
for p, ID in zip(particles,ids):
    p.id = ID
    p.x = normal(loc=0, scale=1.5, size=num_timesteps)
    p.y = normal(loc=0, scale=1.5, size=num_timesteps)
    p.z = normal(loc=0, scale=1.5, size=num_timesteps)

    p.vx = normal(loc=0, scale=0.5, size=num_timesteps)
    p.vy = normal(loc=0, scale=0.5, size=num_timesteps)
    p.vz = normal(loc=0, scale=0.5, size=num_timesteps)

base_p = "halos_snapshot_%i.hdf5.%s.particles"
p_files = [base_p%(i,"%i") for i in range(num_timesteps)]

for pf, t in zip(p_files, range(num_timesteps)):

    for j in range(num_pfiles):
        f = open(pf%j,'w')

        l = p_preamble%(times_list[t])
        f.write(l)

        for p,i in zip(particles,range(len(particles))):

            if (i+j)%num_pfiles == 0:
                line = [p.x[t], p.y[t], p.z[t], p.vx[t], p.vy[t], p.vz[t],\
		        p.id, 0, 0, 0]
                line = [str(l) for l in line]
                f.write(" ".join(line)+'\n')

        f.close()


####################
# write times list #
####################

f = open("output_snapshot_times.txt",'w')

for t in times_list:
    f.write("%s\n"%str(t))



######################################
## write consistent-trees test data ##
######################################

preamble = "#ID DescID Mvir Rvir X Y Z VX VY VZ PID UPID Original_ID\n\
#Units: Masses in Msun / h\n\
#Units: Positions in Mpc / h (comoving)\n\
#Units: Velocities in km / s (physical)\n\
#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n\
#Units: Radii in kpc / h (comoving)\n\
#Units: Tidal force in dimensionless units (Rhalo / Rhill)\n"

ID = [0]*num_timesteps
Mvir = [1233453 for i in range(num_timesteps)]
Rvir = [1231 for i in range(num_timesteps)]
X = normal(loc=0, scale=1, size=num_timesteps)
Y = normal(loc=0, scale=1, size=num_timesteps)
Z = normal(loc=0, scale=1, size=num_timesteps)
VX = normal(loc=0, scale=1, size=num_timesteps)
VY  = normal(loc=0, scale=1, size=num_timesteps)
VZ = normal(loc=0, scale=1, size=num_timesteps)
PID = [0 for i in range(num_timesteps)]
UPID = [0 for _ in range(num_timesteps)]
OID = [3 for _ in range(num_timesteps)]

ID_2 = [1]*num_timesteps
OID_2 = [4 for _ in range(num_timesteps)]

base_rs = "really_consistent_%i.list"
rs_files = [base_rs%i for i in range(num_timesteps)]

for fn, i in zip(rs_files,range(num_timesteps)):

    f = open(fn,'w')

    f.write(preamble)

    line = [ID[i], ID[i], Mvir[i], Rvir[i], \
            X[i], Y[i], Z[i], VX[i], VY[i], VZ[i], PID[i], UPID[i], OID[i]]
    line = [str(l) for l in line]

    f.write(" ".join(line)+'\n')

    line = [ID_2[i], ID_2[i], Mvir[i], Rvir[i], \
            X[i], Y[i], Z[i], VX[i], VY[i], VZ[i], ID[i], UPID[i], OID_2[i]]
    line = [str(l) for l in line]

    f.write(" ".join(line))

    f.close()













