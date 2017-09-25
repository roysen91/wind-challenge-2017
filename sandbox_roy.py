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




# Arrays for ploting
R_L_array = np.array([2,3,10,15,20,25,30]) # resistance load range will be between 0.5 and 50 Ohm (see Mail)
rps_array = np.linspace(1,50,50)

# first wire diameter [mm]
d_kabel = 1.2*10**-3 
rho_cu 	= 15.11 * 10**-9  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf
#Calculation parameters
l_wire 	= fcn.l_kabel(N,l_coil_turn)
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= num_coils * fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
print(B_peak)
L_coil 	= num_coils *fcn.L(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_eff,i_eff,P_prox,P_skin,P_wire = (np.zeros([R_L_array.size,rps_array.size]) for i in range(7))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_eff[j][i]		= fcn.e_eff(B_peak,l_eff_generator,v[j][i],N)
		i_eff[j][i] 	= fcn.i_eff(e_eff[j][i], R_i, R_L, L_coil, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_eff[j][i],d_kabel,f[j][i])*np.sqrt(2) # peak value of amperage
		P_wire[j][i] 	= fcn.i_eff(e_eff[j][i],R_i,R_L,L_coil,f[j][i])

#total loss power and load power
P_loss 	= P_wire#+P_prox+P_skin
P_L 	= fcn.P_L(R_L,i_eff)
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
d_kabel = 1.5*10**-3
rho_cu 	= 9.756 * 10**-9  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf
#Calculation parameters
l_wire 	= fcn.l_kabel(N,l_coil_turn )
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= num_coils *fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
print(B_peak)
L_coil 	= num_coils *fcn.L(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_eff,i_eff,P_prox,P_skin,P_wire = (np.zeros([R_L_array.size,rps_array.size]) for i in range(7))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_eff[j][i]		= fcn.e_eff(B_peak,l_eff_generator,v[j][i],N)
		i_eff[j][i] 	= fcn.i_eff(e_eff[j][i], R_i, R_L, L_coil, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_eff[j][i],d_kabel,f[j][i])*np.sqrt(2) # peak value of amperage
		P_wire[j][i] 	= fcn.i_eff(e_eff[j][i],R_i,R_L,L_coil,f[j][i])

#total loss power and load power
P_loss 	= P_wire#+P_prox+P_skin
P_L 	= fcn.P_L(R_L,i_eff)
M_G		= (P_L+P_loss)/rps_array

M_G = M_G.T

ax2.plot(rps_array,M_G)
ax2.plot(rps_T_N,M_T_N,'ro')
ax2.set_xlabel('rps [1/s]')
ax2.set_ylabel('M_G [Nm]')
ax2.legend([str(R_L) for R_L in R_L_array])
ax2.set_title('1.5 mm')

plt.show()