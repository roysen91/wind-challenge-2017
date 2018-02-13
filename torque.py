import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv


radial_gen = gn.Generator(design='radial',coils=4)

rps_array = np.linspace(1,50,50)
M_rad,eta_rad,R_L_rad = (np.zeros([rps_array.size]) for i in range(3))

v=10
r=450*10**-3
rho=1.225



Turb_n = np.array([])
Turb_M = np.array([])
with open('M_n.txt','r') as file_handle:
	for line in file_handle:
	    a,b = line.strip().split('\t')
	    Turb_n=np.append(Turb_n,float(a))
	    Turb_M=np.append(Turb_M,float(b))


Turb_n = Turb_n/60
Turb_M = Turb_M

P_wind = 1/2*rho*v**3*np.pi*r**2
P_rotor =2*np.pi*Turb_n*Turb_M

cp=P_rotor/P_wind

cp=np.interp(rps_array,Turb_n,cp)

print(cp)



radial_gen.coil.set_d_wire(1.25*10**-3,0.07*10**-3)
radial_gen.coil.set_windings(55)

for i,rps in enumerate(rps_array):
	radial_gen.rot_speed = rps
	radial_gen.M_T = np.interp(rps,Turb_n,Turb_M)
	radial_gen.compute()
	M_rad[i] = radial_gen.M
	eta_rad[i] = radial_gen.eff
	if radial_gen.R_L>50:
		R_L_rad[i] = 0
	else:
		R_L_rad[i] = radial_gen.R_L

eta_ges=cp*eta_rad

fig, [ax1,ax2,ax3] = plt.subplots(3)
ax1.plot(rps_array*60,M_rad.T)
ax1.plot(Turb_n*60,Turb_M,'k')

ax2.plot(rps_array*60,R_L_rad.T)
ax3.plot(rps_array*60,eta_rad.T,rps_array*60,eta_ges)


ax1.set_xlabel('rps [1/min]')
ax2.set_xlabel('rps [1/min]')
ax3.set_xlabel('rps [1/min]')
ax1.set_ylabel('Drehmoment [Nm]')
ax2.set_ylabel('Lastwiderstand [Ohm]')
ax3.set_ylabel('Wirkungsgrad [1]')
ax1.legend(loc="upper right")
ax3.legend(['eta_gen','eta_gesamt'])
plt.show()