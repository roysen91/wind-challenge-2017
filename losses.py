import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn
from axial_gen import *
#from radial_gen import *


# first wire diameter [mm]
R_L = 20

# Arrays for ploting
N = np.array([1,2,3,4,5,6,7,8])
wire = 0


#Calculation parameters
x_mag 	= fcn.x_mag(N,d_kabel[wire],lagen)
B_peak 	= fcn.B_peak(x_mag)
f		= fcn.f(rps_T_N,num_pole_pairs)
v 		= fcn.v_radial(rps_T_N,r_mag)

l_wire 	= fcn.l_kabel(N,l_coil_turn)
A_wire  = fcn.A_kabel(d_kabel[wire])

R_i 	= num_coils * fcn.R_i(rho_cu[wire],l_wire,A_wire)
L_coil 	= num_coils * fcn.Inductance(N, r_outer, r_inner,l_coil_eff,num_pole_pairs)
e_peak	= num_coils * fcn.Voltage(B_peak,l_coil_eff,v,N)
i_peak 	= fcn.Current(e_peak, R_i, R_L, L_coil, f)

# calculate power of losses
P_prox  	= fcn.P_prox(i_peak,d_kabel[wire],f)
P_skin  	= fcn.P_skin(i_peak,d_kabel[wire],f)
P_wire 		= R_i*i_peak**2

fig,ax2 = plt.subplots()

ax2.plot(N,P_wire,N,P_wire+P_skin,N,P_wire+P_skin+P_prox)

ax2.set_xlabel('N [1]')
ax2.set_ylabel('P [W]')
ax2.legend(['P_wire','+P_skin','+P_prox'])

plt.show()