import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

radial_gen_2 = gn.Generator(design='radial',coils=2)
radial_gen_4 = gn.Generator(design='radial',coils=4)

N_array 	= np.linspace(1,200,100)
#N_array 	= np.array([100])
d_array		= np.array([1.25])*10**-3#[mm]
#d_array 	= np.array([0.001])
#RL_array 	= np.linspace(1,50,10)


eta_rad_2,eta_rad_4,R_L_rad_2,R_L_rad_4,l_wire_rad_2,l_wire_rad_4 = (np.zeros([N_array.size,d_array.size]) for i in range(6))

for j,d in enumerate(d_array):
	radial_gen_2.coil.set_d_wire(d,0.07*10**-3)
	radial_gen_4.coil.set_d_wire(d,0.07*10**-3)
	for i,N in enumerate(N_array):
		
		#axial_gen.coil.set_windings(N)
		#axial_gen.compute()
		#eta_ax[i][j] = axial_gen.eff
		#R_L_ax[i][j] = axial_gen.R_L
		#l_wire_ax[i][j] = axial_gen.l_wire

		radial_gen_2.coil.set_windings(N)
		radial_gen_2.compute()
		eta_rad_2[i][j] = radial_gen_2.eff
		R_L_rad_2[i][j] = radial_gen_2.R_L
		l_wire_rad_2[i][j] = radial_gen_2.l_wire


		radial_gen_4.coil.set_windings(N)
		radial_gen_4.compute()
		eta_rad_4[i][j] = radial_gen_4.eff
		R_L_rad_4[i][j] = radial_gen_4.R_L
		l_wire_rad_4[i][j] = radial_gen_4.l_wire
	#print(d,N,radial_gen.coil.layer_x_max,radial_gen.coil.layer_x,radial_gen.coil.layer_y)



fig, (ax1,ax2) = plt.subplots(2)
#ax1.plot(N_array,eta_ax)#,N_array,eta_rad.T)
#ax1.set_xlabel('N [1]')
#ax1.set_ylabel('eta Generator [1]')
#ax1.legend([str(d) for d in d_array*10**3])
#ax1.set_ylim([0,1])
#ax1.grid()

#N = 64.0
#i_d = int(np.where( d_array==1.5*10**-3)[0])
#eta_N = np.round(np.interp(N,N_array,eta_ax[:,i_d]),3)
#ax1.scatter(N, eta_N)
#R_L = np.round(np.interp(N,N_array,R_L_ax[:,i_d]),2)
#txt1 = 'R_L = '+str(R_L)+' Ohm \n'+'eta = '+str(eta_N)
#ax1.annotate(txt1, xy=(N,eta_N))

ax1.plot(N_array,eta_rad_2)
ax1.set_xlabel('N [1]')
ax1.set_ylabel('eta Generator [1]')
ax1.legend([str(d) for d in d_array*10**3])
ax1.set_ylim([0,1])
ax1.grid()

N = 110
i_d = int(np.where( d_array==1.25*10**-3)[0])
eta_N = np.round(np.interp(N,N_array,eta_rad_2[:,i_d]),3)
ax1.scatter(N, eta_N)
R_L = np.round(np.interp(N,N_array,R_L_rad_2[:,i_d]),2)
txt1 = 'R_L = '+str(R_L)+' Ohm \n'+'eta = '+str(eta_N)
ax1.annotate(txt1, xy=(N,eta_N))

ax2.plot(N_array,eta_rad_4)
ax2.set_xlabel('N [1]')
ax2.set_ylabel('eta Generator [1]')
ax2.legend([str(d) for d in d_array*10**3])
ax2.set_ylim([0,1])
ax2.grid()

N = 55
i_d = int(np.where( d_array==1.25*10**-3)[0])
eta_N = np.round(np.interp(N,N_array,eta_rad_4[:,i_d]),3)
ax2.scatter(N, eta_N)
R_L = np.round(np.interp(N,N_array,R_L_rad_4[:,i_d]),2)
txt1 = 'R_L = '+str(R_L)+' Ohm \n'+'eta = '+str(eta_N)
ax2.annotate(txt1, xy=(N,eta_N))


plt.show()