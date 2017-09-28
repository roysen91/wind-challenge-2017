import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

axial_gen = gn.Generator()
radial_gen = gn.Generator(design='radial')

R_L_array = np.array([1,3,5,10,15,20,25,30,40,50]) # resistance load range will be between 0.5 and 50 Ohm (see Mail)
N_array 	= np.array([1,2,3,4,5,6,7,8,9,10])

rps_array = np.linspace(1,50,50)

M_rad,M_ax = (np.zeros([R_L_array.size,rps_array.size]) for i in range(2))

for j,Res in enumerate(R_L_array):
	axial_gen.R_L = Res
	radial_gen.R_L = Res
	for i,rps in enumerate(rps_array):
		
		axial_gen.rot_speed = rps
		axial_gen.compute()
		M_ax[j][i] = axial_gen.M

		radial_gen.rot_speed = rps
		radial_gen.compute()
		M_rad[j][i] = radial_gen.M

fig, ax = plt.subplots()
r_Turb  = 0.45
v_w 	=10
Turb_n  = np.array([i for i in range(12)])*v_w/(np.pi*2*r_Turb)
Turb_M  = np.array([0.5,0.55,0.7,0.85,0.95,1.1,1.1,1.1,0.9,0.6,0.35,0])

ax.plot(rps_array,M_ax.T)#,rps_array,M_rad.T)
ax.plot(Turb_n,Turb_M,'k')
ax.set_xlabel('rps [1/s]')
ax.set_ylabel('Torque Generator [Nm]')
#ax.legend(['axial','radial'])
ax.legend([str(R) for R in R_L_array])
plt.show()