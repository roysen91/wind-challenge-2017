import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

axial_gen = gn.Generator()
radial_gen = gn.Generator(design='radial')

rps_array = np.linspace(1,50,50)
M_rad,M_ax,R_L_ax = (np.zeros([rps_array.size]) for i in range(3))

# turbine characteristic
turb_nom_rps = 24.76
turb_nom_moment = 1.1
r_Turb  = 0.45
v_w 	= 10
Turb_n  = np.array([i for i in range(12)])*v_w/(np.pi*2*r_Turb)
Turb_M  = np.array([0.5,0.55,0.7,0.85,0.95,1.1,1.1,1.1,0.9,0.6,0.35,0])


axial_gen.coil.set_d_wire(1.5*10**-3)
radial_gen.coil.set_d_wire(1.5*10**-3)

axial_gen.coil.set_windings(60)
radial_gen.coil.set_windings(60)

for i,rps in enumerate(rps_array):
	axial_gen.rot_speed = rps
	axial_gen.M_T = np.interp(rps,Turb_n,Turb_M)
	axial_gen.compute()
	M_ax[i] = axial_gen.M
	R_L_ax[i]= axial_gen.R_L

	radial_gen.rot_speed = rps
	radial_gen.compute()
	M_rad[i] = radial_gen.M

fig, ax = plt.subplots()
ax.plot(rps_array,M_ax.T)#,rps_array,M_rad.T)
ax.plot(Turb_n,Turb_M,'k')


ax.scatter(turb_nom_rps, turb_nom_moment)
txt = 'R_L = '+str(np.round(np.interp(turb_nom_rps,rps_array,M_ax),2))+' Ohm'
ax.annotate(txt, (turb_nom_rps,turb_nom_moment))

ax.set_xlabel('rps [1/s]')
ax.set_ylabel('Torque Generator [Nm]')
#ax.legend(['axial','radial'])
#ax.legend([str(R) for R in R_L_ax])
print(R_L_ax)
plt.show()