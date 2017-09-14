import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn



# Stellgrößen (Werden nur verwendet, wenn diese konstant gehalten werden sollen):
num_pole_pairs 		= 2  				# Polpaarzahl
N 					= 20  # Wicklungen
lagen 				= 10  # Zahl der Kabelschichten
d_kabel 			= 0.00175  # Kabeldurchmesser
rho_cu 				= 17 * 10**-9  # sp. Widerstand Kupfer [Ohm+m]


# Variablen
M_T 				= 6.924  # Moment Welle [Nm]
r_mag 				= 67.5 * 10**-3  # Radius Wellenmitte bis Magnetmitte
l_eff 				= 2 * 0.045	
l_coil_aus 			= 0.9 * 2 * np.pi * 120 / 360 / num_pole_pairs  # Außendurchmesser Spule
l_coil_inn 			= 0.45 * 2 * np.pi * 120 / 360 / num_pole_pairs  # Innendurchmesser der Spule
P_w 				= 1.2 / 2 * 0.45**2 * np.pi * 10**3  # dichte luft = 1.2, rotorradius = 0.45, windgeschwindigkeit = 10




R_L_array 				= np.array([1,3,5,10,15,20,25,30])
#rps 				= 1485 / 60  # Umdrehungen d. Welle [1/s]
rps_array = np.linspace(1,100,50)

# first wire diameter [mm]
d_kabel = 1.2*10**-3 
#Calculation parameters
l_wire 	= fcn.l_kabel(N,num_pole_pairs,l_coil_inn, l_coil_aus, l_eff)
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
L 		= fcn.L(N,l_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_eff,i_eff,P_prox,P_skin,P_wire = (np.zeros([R_L_array.size,rps_array.size]) for i in range(7))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_eff[j][i]		= fcn.e_eff(B_peak,l_eff,v[j][i],N)
		i_eff[j][i] 	= fcn.i_eff(e_eff[j][i], R_i, R_L, L, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_eff[j][i],d_kabel,f[j][i])*np.sqrt(2) # peak value of amperage
		P_wire[j][i] 	= fcn.i_eff(e_eff[j][i],R_i,R_L,L,f[j][i])

#total loss power and load power
P_loss 	= P_wire#+P_prox+P_skin
P_L 	= fcn.P_L(R_L,i_eff)
M_T		= (P_L+P_loss)/rps_array

M_T = M_T.T

#plot torque-rpm diagram
fig, (ax1,ax2) = plt.subplots(2,sharey=True)
ax1.plot(rps_array,M_T)
ax1.set_xlabel('rps [1/s]')
ax1.set_ylabel('M_T [Nm]')
ax1.legend([str(R_L) for R_L in R_L_array])
ax1.set_title('1.2 mm')

# second wire diameter [mm]
d_kabel = 1.5*10**-3 
#Calculation parameters
l_wire 	= fcn.l_kabel(N,num_pole_pairs,l_coil_inn, l_coil_aus, l_eff)
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B_peak 	= fcn.B_peak(x_mag)
L 		= fcn.L(N,l_eff,num_pole_pairs)

# initiate np.arrays 
f,v,e_eff,i_eff,P_prox,P_skin,P_wire = (np.zeros([R_L_array.size,rps_array.size]) for i in range(7))

for j,R_L in enumerate(R_L_array):
	for i,rps in enumerate(rps_array):
		f[j][i] 		= fcn.f(rps,num_pole_pairs)
		v[j][i] 		= fcn.v_radial(rps,r_mag)
		e_eff[j][i]		= fcn.e_eff(B_peak,l_eff,v[j][i],N)
		i_eff[j][i] 	= fcn.i_eff(e_eff[j][i], R_i, R_L, L, f[j][i])


		# calculate power of losses
		P_prox[j][i]  	= fcn.P_prox(d_kabel,f[j][i])
		P_skin[j][i]  	= fcn.P_skin(i_eff[j][i],d_kabel,f[j][i])*np.sqrt(2) # peak value of amperage
		P_wire[j][i] 	= fcn.i_eff(e_eff[j][i],R_i,R_L,L,f[j][i])

#total loss power and load power
P_loss 	= P_wire#+P_prox+P_skin
P_L 	= fcn.P_L(R_L,i_eff)
M_T		= (P_L+P_loss)/rps_array

M_T = M_T.T

ax2.plot(rps_array,M_T)
ax2.set_xlabel('rps [1/s]')
ax2.set_ylabel('M_T [Nm]')
ax2.legend([str(R_L) for R_L in R_L_array])
ax2.set_title('1.5 mm')

plt.show()