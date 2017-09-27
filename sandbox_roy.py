import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn
from axial_gen import *
#from radial_gen import *

# Arrays for ploting
R_L_array = np.array([1,3,5,10,15,20,25,30]) # resistance load range will be between 0.5 and 50 Ohm (see Mail)
rps_array = np.linspace(1,50,50)

# first wire diameter [m]
d_kabel = 1.2
rho_cu 	= 15.11 * 10**-3  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf
#Calculation parameters
l_wire 	= fcn.l_kabel(N,l_coil_turn)
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= num_coils * fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
L_coil 	= num_coils *fcn.Inductance(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_peak,i_peak,P_prox,P_skin,P_wire,P_L = (np.zeros([R_L_array.size,rps_array.size]) for i in range(8))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_peak[j][i]	= num_coils * fcn.Voltage(B_peak,l_coil_eff,v[j][i],N)
		i_peak[j][i] 	= fcn.Current(e_peak[j][i], R_i, R_L, L_coil, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(i_peak[j][i],d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_peak[j][i],d_kabel,f[j][i])
		P_wire[j][i] 	= R_i*i_peak[j][i]**2

		P_L[j][i] 	= fcn.P_L(R_L,i_peak[j][i])

#total loss power and load power
P_loss 	= P_prox+P_skin
M_G		= (P_L+P_loss)/rps_array

M_G = M_G.T

#plot torque-rpm diagram
fig, (ax1,ax2) = plt.subplots(2,sharey=True)
ax1.plot(rps_array,M_G)
ax1.plot(rps_T_N,M_T_N,'ro')
ax1.set_xlabel('rps [1/s]')
ax1.set_ylabel('M_G [Nm]')
ax1.legend([str(R_L) for R_L in R_L_array])
ax1.set_title('1.2 mm')

# second wire diameter [mm]
d_kabel = 1.5
rho_cu 	= 9.756 * 10**-3  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf
#Calculation parameters
l_wire 	= fcn.l_kabel(N,l_coil_turn )
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= num_coils *fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
L_coil 	= num_coils *fcn.Inductance(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_peak,i_peak,P_prox,P_skin,P_wire,P_L = (np.zeros([R_L_array.size,rps_array.size]) for i in range(8))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_peak[j][i]	= num_coils * fcn.Voltage(B_peak,l_coil_eff,v[j][i],N)
		i_peak[j][i] 	= fcn.Current(e_peak[j][i], R_i, R_L, L_coil, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(i_peak[j][i],d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_peak[j][i],d_kabel,f[j][i])
		P_wire[j][i] 	= R_i*i_peak[j][i]

		P_L[j][i] 	= fcn.P_L(R_L,i_peak[j][i])

#total loss power and load power
P_loss 	= P_prox+P_skin
M_G		= (P_L+P_loss)/rps_array

M_G = M_G.T

ax2.plot(rps_array,M_G)
ax2.plot(rps_T_N,M_T_N,'ro')
ax2.set_xlabel('rps [1/s]')
ax2.set_ylabel('M_G [Nm]')
ax2.legend([str(R_L) for R_L in R_L_array])
ax2.set_title('1.5 mm')

plt.show()