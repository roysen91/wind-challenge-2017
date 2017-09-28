import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

axial_gen = gn.Generator()
radial_gen = gn.Generator(design='radial')

N_array 	= np.array([1,2,3,4,5,6,7,8,9,10])
d_array		= np.array([1.2, 1.5 ]) #[mm]
rho_array 	= np.array([15.11 * 10**-3, 9.756 * 10**-3]) 


P_wire_ax,P_wire_rad,P_vortex_ax,P_vortex_rad = (np.zeros([d_array.size,N_array.size]) for i in range(4))

for j,d in enumerate(d_array):
	axial_gen.coil.d_wire = d
	axial_gen.coil.rho_mat = rho_array[j]
	radial_gen.coil.d_wire = d
	radial_gen.coil.rho_mat = rho_array[j]
	
	for i,N in enumerate(N_array):
		
		axial_gen.coil.windings = N
		axial_gen.compute()
		P_wire_ax[j][i] = axial_gen.P_wire
		P_vortex_ax[j][i] = axial_gen.P_skin+axial_gen.P_prox

		radial_gen.coil.windings = N
		radial_gen.compute()
		P_wire_rad[j][i] = radial_gen.P_wire
		P_vortex_rad[j][i] = radial_gen.P_skin+radial_gen.P_prox

print()
fig, ax = plt.subplots()
ax.plot(N_array,P_wire_ax.T+P_vortex_ax.T,N_array,P_wire_rad.T+P_vortex_rad.T)
ax.set_xlabel('N [1]')
ax.set_ylabel('P_loss [W]')
ax.legend(['axial 1.2 mm','axial 1.5 mm','radial 1.2 mm','radial 1.5 mm'])
plt.show()