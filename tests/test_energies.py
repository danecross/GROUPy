
from Particle import *
from Energies import *
from Halo import *

m = 4.61528e+06

# create test "field"
p_field = [Particle(0,0,0,0,0,0,mass=m), Particle(1,0,0,1,0,0,mass=m)]

h = Halo()
h.vx = h.vy = h.vz = [0]

##############
## test PEg ##
##############

PE = PEg(p_field[0], p_field, 0)
assert(PE == -9.618380519971521e-35)

#############
## test KE ##
#############

KE = KE(p_field[1],0)
assert(KE == 2.4231272121835527e-33)



