import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn

# Constants
num_pole_pairs 		= 2  				# Polpaarzahl
num_coils = 4 # number of coils used [1]

# Variable generator parameter
N 					= 20  # Wicklungen
lagen 				= 10  # Zahl der Kabelschichten
d_kabel 			= 0.00175  # Kabeldurchmesser
# Radialgenerator
r_inner 			= 45.5*10**-3 # inner radius of magnets [mm]
r_outer 			= 90.5*10**-3 # outer radius of magnets [mm]
r_mag 				= (r_inner+r_outer)/2  # Radius Wellenmitte bis Magnetmitte [m]

l_coil_eff 			= 0.045		# effective coil length per turn [m]
l_coil_outer 		= (r_outer * 2 * np.pi / 360) * 60  # outer coil length per turn [m]
l_coil_inner 		= (r_inner * 2 * np.pi / 360) * 60# inner coil length per turn [m]
l_coil_turn 		= 2*l_coil_eff + l_coil_outer + l_coil_inner

l_eff_generator = 2 * l_coil_eff * num_coils


# Turbine data
M_T_N 	= 1.11 # nominal turbine torque [Nm]
rps_T_N = 24.76 #nominal turbine rpm [1/s]
R_L = 10 # load resistance [Ohm]


# Arrays for ploting
R_L_array = np.array([2,3,10,15,20,25,30]) # resistance load range will be between 0.5 and 50 Ohm (see Mail)
rps_array = np.linspace(1,50,50)
N_array = np.linspace(20,50,50)

# first wire diameter [mm]
d_kabel_array = np.array([1.2*10**-3, 1.5*10**-3 ])
rho_cu 	= np.array([15.11 * 10**-9, 9.756 * 10**-9])  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf

#Calculation parameters


# initiate np.arrays 
l_wire,A_wire,R_i,x_mag,B_peak,L_coil,f,v,e_eff,i_eff,P_prox,P_skin,P_wire,eta_G,P_loss = (np.zeros([d_kabel_array.size,N_array.size]) for i in range(15))

for j,d_kabel in enumerate(d_kabel_array):
	for i,N in enumerate(N_array):

		l_wire[j][i] 	= fcn.l_kabel(N,l_coil_turn)
		A_wire[j][i]  	= fcn.A_kabel(d_kabel)
		R_i[j][i] 		= num_coils * fcn.R_i(rho_cu[j],l_wire[j][i],A_wire[j][i])
		x_mag[j][i] 	= fcn.x_mag(N,d_kabel,lagen)
		B_peak[j][i] 	= fcn.B_peak(x_mag[j][i])
		L_coil[j][i] 	= num_coils *fcn.L(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)

		f[j][i] 		= fcn.f(rps_T_N,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps_T_N,r_mag)
		e_eff[j][i]		= fcn.e_eff(B_peak[j][i],l_eff_generator,v[j][i],N)
		i_eff[j][i] 	= fcn.i_eff(e_eff[j][i], R_i[j][i], R_L, L_coil[j][i], f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_eff[j][i],d_kabel,f[j][i])*np.sqrt(2) # peak value of amperage
		P_wire[j][i] 	= fcn.i_eff(e_eff[j][i],R_i[j][i],R_L,L_coil[j][i],f[j][i])

		#total loss power and load power
		P_loss[j][i] 	= P_wire[j][i]#+P_prox+P_skin
		eta_G[j][i]		=1/(1+P_loss[j][i]/(e_eff[j][i]*i_eff[j][i]))

eta_G = eta_G.T

#plot torque-rpm diagram
fig, ax = plt.subplots()
ax.plot(N_array,eta_G)
ax.set_xlabel('N [1]')
ax.set_ylabel('eta Generator [1]')
ax.legend([str(d) for d in d_kabel_array])
plt.show()
