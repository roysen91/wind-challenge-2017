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




#R_L_array 				= np.array([1,5,10,20,30,40,50])
R_L = 10
#rps 				= 1485 / 60  # Umdrehungen d. Welle [1/s]
rps_array = np.linspace(1,2000,50)



#Calculation of losses
l_wire 	= fcn.l_kabel(N,num_pole_pairs,l_coil_inn, l_coil_aus, l_eff)
A_wire  = fcn.A_kabel(d_kabel)
R_i 	= fcn.R_i(rho_cu,l_wire,A_wire)
x_mag 	= fcn.x_mag(N,d_kabel,lagen)
B 		= fcn.B_field(x_mag)
L 		= fcn.L(N,l_eff,num_pole_pairs)


f,v,e_indu,i_indu,P_prox,P_skin,P_wire = (np.zeros([rps_array.size]) for i in range(7))

for i,rps in enumerate(rps_array):
	f[i] 		= fcn.f(rps,num_pole_pairs)
	v[i] 		= fcn.v_radial(rps,r_mag)
	e_indu[i]	= fcn.e_indu(B,l_eff,v[i],N)
	i_indu[i] 	= fcn.i_indu(e_indu[i], R_i, R_L, L, f[i])


	# P_loss = I_peak*R_wire + P_Prox + P_skin)
	P_prox[i]  	= fcn.P_prox(d_kabel,f[i])
	P_skin[i]  	= fcn.P_skin(i_indu[i]*np.sqrt(2),d_kabel,f[i])
	P_wire[i] 	= fcn.i_indu(e_indu[i],R_i,R_L,L,f[i])*np.sqrt(2)

P_loss 	= P_prox+P_skin+P_wire
P_L 	= fcn.P_L(R_L,i_indu*np.sqrt(2))
M_T		= (P_L+P_loss)/rps_array


plt.plot(rps_array,M_T)
plt.xlabel('rps [1/s]')
plt.ylabel('M_T [Nm]')
plt.show()