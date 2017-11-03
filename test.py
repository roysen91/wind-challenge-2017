import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


radial_gen_4 = gn.Generator(design='radial',coils=4)
radial_gen_4.coil.set_d_wire(1.25*10**-3,0.07*10**-3)
radial_gen_4.coil.set_windings(55)
radial_gen_4.compute()

radial_gen_2 = gn.Generator(design='radial',coils=2)
radial_gen_2.coil.set_d_wire(1.25*10**-3,0.07*10**-3)
radial_gen_2.coil.set_windings(110)
radial_gen_2.compute()

print(radial_gen_4)
#print(radial_gen_2)
#print(radial_gen_4.coil.max_coil_width)