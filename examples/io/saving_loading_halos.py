
import os
from Halo import *

# this is the path that we will save the halos to
save_dir = "save_output/"
if not os.path.exists(save_dir): os.mkdir(save_dir)
save_name = os.path.join(save_dir, "halo.npy")

# create Halo object

h = Halo()

# add some fake information
h.ids = [1, 2, 3, 4]
h.x = [1, 2, 3, 4]

# save the halo

h.save(save_name)

# re-load the halo

h2 = load_halo(save_name)

# this should print True
print(h2.ids == h.ids)






