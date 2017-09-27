import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn
from axial_gen import *
#from radial_gen import *


# Arrays for ploting
R_L = 20 
N_array = np.array([1,2,3,4,5,6,7,8])

# initiate np.arrays 
l_wire,A_wire,R_i,x_mag,B_peak,L_coil,f,v,e_peak,i_peak,P_prox,P_skin,P_wire,eta_G,P_loss = (np.zeros([d_kabel.size,N_array.size]) for i in range(15))

for j,d_wire in enumerate(d_kabel):
	for i,N in enumerate(N_array):

		l_wire[j][i] 	= fcn.l_kabel(N,l_coil_turn)
		A_wire[j][i]  	= fcn.A_kabel(d_wire)
		x_mag[j][i] 	= fcn.x_mag(N,d_wire,lagen)
		B_peak[j][i] 	= fcn.B_peak(x_mag[j][i])
		f[j][i] 		= fcn.f(rps_T_N,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps_T_N,r_mag)


		R_i[j][i] 		= num_coils * fcn.R_i(rho_cu[j],l_wire[j][i],A_wire[j][i])
		L_coil[j][i] 	= num_coils *fcn.Inductance(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)
		e_peak[j][i]	= num_coils * fcn.Voltage(B_peak[j][i],l_coil_eff,v[j][i],N)
		i_peak[j][i] 	= fcn.Current(e_peak[j][i], R_i[j][i], R_L, L_coil[j][i], f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(i_peak[j][i],d_wire,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_peak[j][i],d_wire,f[j][i])
		P_wire[j][i] 	= R_i[j][i]*i_peak[j][i]**2

		#total loss power and load power
		P_loss[j][i] 	= P_prox[j][i]+P_skin[j][i]+P_wire[j][i]
		eta_G[j][i]		= 1/(1+P_loss[j][i]/(e_peak[j][i]*i_peak[j][i]))

eta_G = eta_G.T


#plot torque-rpm diagram
fig, ax = plt.subplots()
ax.plot(N_array,eta_G)
ax.set_xlabel('N [1]')
ax.set_ylabel('eta Generator [1]')
ax.legend([str(d) for d in d_kabel])
plt.show()
