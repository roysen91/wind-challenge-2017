import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

axial_gen = gn.Generator()
radial_gen = gn.Generator(design='radial')

N_array 	= np.array([1,2,3,4,5,6,7,8,9,10])
d_array		= np.array([1.2, 1.5 ]) #[mm]
rho_array 	= np.array([15.11 * 10**-3, 9.756 * 10**-3]) 


eta_rad,eta_ax = (np.zeros([d_array.size,N_array.size]) for i in range(2))

for j,d in enumerate(d_array):
	axial_gen.coil.d_wire = d
	axial_gen.coil.rho_mat = rho_array[j]
	radial_gen.coil.d_wire = d
	radial_gen.coil.rho_mat = rho_array[j]
	
	for i,N in enumerate(N_array):
		
		axial_gen.coil.windings = N
		axial_gen.compute()
		eta_ax[j][i] = axial_gen.eff

		radial_gen.coil.windings = N
		radial_gen.compute()
		eta_rad[j][i] = radial_gen.eff


print(radial_gen.P_L,radial_gen.V*radial_gen.I,radial_gen.P_loss)

fig, ax = plt.subplots()
ax.plot(N_array,eta_ax.T,N_array,eta_rad.T)
ax.set_xlabel('N [1]')
ax.set_ylabel('eta Generator [1]')
ax.legend(['axial 1.2 mm','axial 1.5 mm','radial 1.2 mm','radial 1.5 mm'])
plt.show()